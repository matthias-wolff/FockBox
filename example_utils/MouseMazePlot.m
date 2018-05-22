classdef MouseMazePlot < handle
  % A 3D mouse-maze plot.
  %
  % preview 3D objecs:
  %   MouseMazePlot.preview(MouseMazePlot.getFloorTileModel())
  %   MouseMazePlot.preview(MouseMazePlot.getWallModel())
  %   MouseMazePlot.preview(MouseMazePlot.getMouseModel())
  %   MouseMazePlot.preview(MouseMazePlot.getCheeseModel())
  % 
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %
  % TODO:
  % - Improve arrangement of multiple scene objects at the same location
  %   -> method arrangeSceneObjects(obj,x,y)
  % - Implement generic scene object handlers addX, removeA, etc.
  
  %% == Constants ==
  
  properties(Constant)
    
    clrFloor   = [0 0.61 0.82];                                                 % Floor edge/lines: KT blue
    %clrAmbient = [0.5*249/255 0.5*129/255 0.5*67/255];                         % Ambient light: orange
    clrAmbient = [0.5 0.5 0.5];                                                 % Ambient light: grey
    cacheDir   = [ mfilename('fullpath') '_cache' ];                            % Cache folder
    
  end
  
  %% == Read-only properties ==
  
  properties(SetAccess=protected)
    
    % The axes of this plot.
    ax;

    % Scene objects.
    sco;

  end
  
  %% == Public API ==
  
  methods

    function obj=MouseMazePlot(arg1)
      % Creates a mouse-maze plot.
      %
      %   obj.MouseMazePlot(mdim)  - Create a mouse-maze plot
      %   MouseMazePlot CLEARCACHE - Clear patch cache
      %
      % arguments:
      %   mdim - A vector specifying the maze dimensions.
      %
      % returns:
      %   obj  - The mouse-maze plot.
      %
      % example:
      %
      %   mmp = MouseMazePlot([4 3]);
      
      % Detect special arguments                                                % -------------------------------------
      if nargin>0&&(ischar(arg1)||isstring(arg1))                               % Argument is symbolic >>
        if strcmpi(arg1,'CLEARCACHE')                                           %   Argument is 'CLEARCACHE' >>
          MouseMazePlot.clearCache();                                           %     Clear patch cache
          return;                                                               %     Exit
        elseif strcmpi(arg1,'DEBUG')                                            %   << Argument is 'DEBUG' >>
          MouseMazePlot.DEBUG();                                                %     Invoke debug method
          return;                                                               %     Exit
        else                                                                    %   << Other sybolic argument >>
          error('Unknown option ''%s''.',arg1);                                 %     Dunno...
        end                                                                     %   <<
      end                                                                       % <<
      
      % Initialize maze and scene objects                                       % -------------------------------------
      if nargin==0; arg1 = [ 4 3 ]; end                                         % No argument -> standard maze
      assert(isvector(arg1)&&isnumeric(arg1)&&length(arg1)==2, ...              % Check argument
        'arg1 must be a numeric vector with two elements');                     % ...
      obj.sco.floor   = gobjects(0,0);                                          % Create floor tile objects array
      obj.sco.xlabels = gobjects(0,0);                                          % Create x labels objects array
      obj.sco.ylabels = gobjects(0,0);                                          % Create y labels objects array
      obj.sco.walls   = gobjects(0,0,0);                                        % Create east walls objects array
      obj.sco.mice    = gobjects(0,0);                                          % Create mouse objects array
      obj.sco.cats    = gobjects(0,0);                                          % Create cat objects array
      obj.sco.cheese  = gobjects(0,0);                                          % Create cheese objects array
      obj.sco.water   = gobjects(0,0);                                          % Create water objects array
      
      % Create scene                                                            % -------------------------------------
      obj.ax = gca;                                                             % Remember axes
      axis equal tight off;                                                     % Parametrize plot axes
      set(gcf,'color','w');                                                     % Set figure plot background to white
      for x=1:arg1(1)                                                           % Loop over x positions >>
        for y=1:arg1(2)                                                         %   Loop over y positions >>
          obj.addFloorTile(x,y);                                                %     Add a floor tile
        end                                                                     %   <<
      end                                                                       % <<
      obj.updateLabels();                                                       % Update base labels
      
      % Lighting and camera                                                     % -------------------------------------
      camproj(obj.ax,'perspective');                                            % Set camera to perspective mode
      campos(obj.ax,[-30 20 -10]);                                              % Set camera location
      camup(obj.ax,[0 1 0]);                                                    % Set camera's upright axis
      obj.sco.ambiLight = camlight(obj.ax,'right');                             % Create directed camera ambient light
      obj.sco.ambiLight.Color = MouseMazePlot.clrAmbient;                       % ...
      obj.sco.ambiLight.Style = 'infinite';                                     % ...
      obj.sco.camLight = camlight(obj.ax,'headlight');                          % Create white camera headlight
      obj.sco.camLight.Style = 'infinite';                                      % ...
      cameratoolbar(gcf,'SetCoordSys','y');                                     % Add cam. toolbar, y-axis principal
    end

    function setCamPos(obj,x,y,z)
      % Sets a new camera position
      %
      %    obj.setCamPos(x,y,z)
      %
      % arguments:
      %    x - x coordinate of new camera position.
      %    y - y coordinate of new camera position.
      %    z - new camera elevation.
      campos(obj.ax,[y z x]);
    end
    
    function [xdim,ydim]=getDim(obj,ndim)
      % Returns the dimension of the maze.
      %
      %   [xdim,ydim]=obj.getDim()
      %   dim=obj.getDim(ndim)
      %
      % arguments:
      %   obj  - The mouse-maze plot.
      %   ndim - The direction, 1 for x dimension, 2 for y dimension.
      %
      % returns:
      %   xdim - The x dimension of the maze.
      %   ydim - The y dimension of the maze.
      %   dim  - The dimension in the selected direction.
      if nargin==1
        xdim = size(obj.sco.floor,2);
        ydim = size(obj.sco.floor,1);
      else
        if ndim==1
          xdim = size(obj.sco.floor,2);
        elseif ndim==2
          xdim = size(obj.sco.floor,1);
        else
          error('Argument ''dim'' must be either 1 or 2.');
        end
      end
    end

    function clear(obj)
      % Clears the plot.
      %
      %    obj.clear();
      %
      % arguments:
      %    obj - The mouse-maze plot.
      
      obj.removeAllFloorTiles();
      obj.removeAllWalls();
      obj.removeAllMice();
      obj.removeAllCats();
      obj.removeAllCheese();
      obj.removeAllWater();
    end
    
    % - Floor tiles

    function a=addFloorTile(obj,x,y,p,color)
      % Places a floor tile at a location in the maze. If there is already
      % a floor tile at the specified location, the method will replace it.
      % If the specified location is greater than the current maze limits in 
      % either direction, the method will enlarge the maze.
      % 
      %   a=obj.addFloorTile(x,y)
      %   a=obj.addFloorTile(x,y,p)
      %   a=obj.addFloorTile(x,y,p,color)
      %
      % arguments:
      %   obj   - The mouse-maze plot.
      %   x     - x coordinate of maze location.
      %   y     - y coordinate of maze location.
      %   p     - Probability of floor tile, 0 to remove, default is 1.
      %   color - Color of floor tile face as RGB triplet, default is white.
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      assert(x>0&&y>0,'Maze coordinates must be positive');                     % Check x and y coordinates
      if nargin<4; p = 1; end                                                   % Default probabilty: 1
      if nargin<5                                                               % No color specified >> 
        color = MouseMazePlot.clrFloor+0.75*(1-MouseMazePlot.clrFloor);         %   Def. color: very light KT blue
      end                                                                       % <<
      if obj.isFloorTileAt(x,y)                                                 % There is already a floor tile >>
        delete(obj.sco.floor(y,x));                                             %   Delete from scene
        obj.sco.floor(y,x)=gobjects(1,1);                                       %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove tile
      if p>1; p=1; end                                                          % Limit probability to 1
      c = MouseMazePlot.getFloorTileModel();                                    % Get copy floor tile model
      c.fvcd(1:4,:) = repmat(color,4,1);                                        % Set tile face color
      c.v = c.v+[y 0 x];                                                        % Translate to maze coordinates
      a = patch(obj.ax,'vertices',c.v,'faces',c.f.v,'FaceVertexCData',c.fvcd);  % Create new cheese scene object
      f = fields(c.p); for i=1:length(f); a.set(f{i},c.p.(f{i})); end           % Copy patch properties
      a.FaceAlpha = p;                                                          % Probability -> face alpha
      a.EdgeAlpha = p;                                                          % Probability -> edge alpha
      obj.sco.floor(y,x) = a;                                                   % Store scene object
      obj.updateLabels();                                                       % Update row and column labels
      camva(obj.ax,'auto');                                                     % Update camera viewing angle
    end

    function removeFloorTile(obj,x,y)
      % Removes a floor tile. If there is no floor tile at the specified
      % location, the method will do nothing.
      %
      %   obj.removeFloorTile(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addFloorTile(x,y,0);                                                  % Invoke addFloorTile with p=0
    end

    function p=isFloorTileAt(obj,x,y)
      % Tests if there is a floor tile at the specified location.
      % If there is no floor tile or if the specified location is outside the 
      % maze boudaries, the method will return 0.
      %
      %   p = obj.isFloorTileAt(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % returns:
      %   p   - The probability of a floor tile being at the specified location.
      p = 0;                                                                    % Default return value: 0
      if x<=0 || y<=0 || x>size(obj.sco.floor,2) || y>size(obj.sco.floor,1)     % Coordinates outside maze >>
        return;                                                                 %   No tile at (x,y)
      end                                                                       % <<
      if isa(obj.sco.floor(y,x),'matlab.graphics.primitive.Patch')              % Floor tile patch at (x,y) >>
        p = obj.sco.floor(y,x).FaceAlpha;                                       %   Probability is alpha value
      end                                                                       % <<
    end

    function f=getAllFloorTiles(obj)
      % Gets all floor tiles.
      %
      %   f=obj.getAllFloorTiles()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   f   - A struct vector containing the locations f(i).x, f(i).y, the
      %         probabilities f(i).p, and the colors f(i).color of all floor
      %         tiles in the maze.  
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = xdim*ydim;                                                            % Maximal number of floor tiles
      f = struct('x', cell(1,n), 'y', cell(1,n), 'p', cell(1,n), ...            % Pre-allocate result
        'color', cell(1,n));                                                    % ...
      i = 1;                                                                    % Pointer into result
      for x=1:xdim                                                              % Loop over x coordinates >>
        for y=1:ydim                                                            %   Loop over y coordinates >>
          p = obj.isFloorTileAt(x,y);                                           %     Get floor tile prob. at location
          if p                                                                  %     If there is a floor tile >>
            f(i).x = x;                                                         %       Store x coordinate in result
            f(i).y = y;                                                         %       Store y coordinate in result
            f(i).p = p;                                                         %       Store probability in result
            f(i).color = obj.sco.floor(y,x).FaceVertexCData(1,:);               %       Store floor tile color in res.
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      f = f(1:i-1);                                                             % Shrink result
    end

    function removeAllFloorTiles(obj)
      % Removes all floor tiles.
      % 
      %   obj.removeAllFloorTiles()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      for i=1:size(obj.sco.floor,2)
        for j=1:size(obj.sco.floor,1)
          obj.addFloorTile(i,j,0);
        end
      end
    end

    % - Walls

    function a=addWall(obj,x,y,d,p,color)
      % Places a wall at a location in the maze. If there is already a wall at 
      % the specified location, the method will replace it. If the specified 
      % location is greater than the current maze limits in either direction, 
      % the method will enlarge the maze.
      % 
      %   a=obj.addWall(x,y,d)
      %   a=obj.addWall(x,y,d,p)
      %   a=obj.addWall(x,y,d,p,color)
      %
      % arguments:
      %   obj   - The mouse-maze plot.
      %   x     - x coordinate of maze location.
      %   y     - y coordinate of maze location.
      %   d     - The wall position, one of 'N', 'E', 'S', or 'W'.
      %   p     - Probability of floor tile, 0 to remove, default is 1.
      %   color - Color of wall faces as RGB triplet, default is white.
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      p0 = obj.isWallAt(x,y,d);                                                 % Get current probability of wall
      [x,y,z] = MouseMazePlot.convertWallCoords(x,y,d);                         % Convert coordinates and direction
      if nargin<5; p = 1; end                                                   % Default probabilty: 1
      if nargin<6                                                               % No color specified >> 
        color = MouseMazePlot.clrFloor+0.75*(1-MouseMazePlot.clrFloor);         %   Def. color: very light KT blue
      end                                                                       % <<
      if p0                                                                     % There is already a wall >>
        delete(obj.sco.walls(y,x,z));                                           %   Delete from scene
        obj.sco.walls(y,x,z)=gobjects(1,1);                                     %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove wall
      if p>1; p=1; end                                                          % Limit probability to 1
      c = MouseMazePlot.getWallModel();                                         % Get copy wall model
      c.fvcd(1:8,:) = repmat(color,8,1);                                        % Set wall face color
      c.v = c.v+[-0.5 0 0];                                                     % Translate to maze coordinates      
      if z==1; c.v = c.v * MouseMazePlot.rotmat(90); end                        % Rotate north and south walls
      c.v = c.v+[y 0 x];                                                        % Translate to maze coordinates
      a = patch(obj.ax,'vertices',c.v,'faces',c.f.v,'FaceVertexCData',c.fvcd);  % Create new cheese scene object
      f = fields(c.p); for i=1:length(f); a.set(f{i},c.p.(f{i})); end           % Copy patch properties
      a.FaceAlpha = p;                                                          % Probability -> face alpha
      a.EdgeAlpha = p;                                                          % Probability -> edge alpha
      obj.sco.walls(y,x,z) = a;                                                 % Store scene object
    end

    function removeWall(obj,x,y,d)
      % Removes a wall. If there is no wall at the specified location, the 
      % method will do nothing.
      %
      %   obj.removeWall(x,y,d)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   d   - The wall position, one of 'N', 'E', 'S', or 'W'.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addWall(x,y,d,0);                                                     % Invoke addWall with p=0
    end

    function p=isWallAt(obj,x,y,d)
      % Tests if there is a wall at the specified location. If there is no wall
      % or if the specified location is outside the maze boudaries, the method 
      % will return 0.
      %
      %   p = obj.isWallAt(x,y,d)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   d   - The wall position, one of 'N', 'E', 'S', or 'W'.
      %
      % returns:
      %   p   - The probability of a floor tile being at the specified location.
      p = 0;                                                                    % Default return value: 0
      [x,y,z] = MouseMazePlot.convertWallCoords(x,y,d);                         % Convert coordinates and direction
      if x<=0 || y<=0 || x>size(obj.sco.walls,2) || y>size(obj.sco.walls,1)     % Coordinates outside maze >>
        return;                                                                 %   No wall at (x,y,d)
      end                                                                       % <<
      if z>size(obj.sco.walls,3); return; end                                   % Check z coordinate
      if isa(obj.sco.walls(y,x,z),'matlab.graphics.primitive.Patch')            % Wall patch at (x,y,d) >>
        p = obj.sco.walls(y,x,z).FaceAlpha;                                     %   Probability is alpha value
      end                                                                       % <<
    end

    function w=getAllWalls(obj)
      % Gets all walls.
      %
      %   w=obj.getAllWalls()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   w   - A struct vector containing the locations w(i).x, w(i).y, the
      %         directions w(i).d, the probabilities w(i).p, and the colors 
      %         w(i).color of all walls in the maze.
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = (xdim+1)*(ydim+1)*2;                                                  % Maximal number of walls
      w = struct('x', cell(1,n), 'y', cell(1,n), 'd', cell(1,n), ...            % Pre-allocate result
        'p', cell(1,n), 'color', cell(1,n));                                    % ...
      i = 1;                                                                    % Pointer into result
      for x=1:xdim+1                                                            % Loop over x coordinates >>
        for y=1:ydim+1                                                          %   Loop over y coordinates >>
          p = obj.isWallAt(x,y,'W');                                            %     Get west wall prob. at location
          if p                                                                  %     If there is a west wall >>
            w(i).x = x;                                                         %       Store x coordinate in result
            w(i).y = y;                                                         %       Store y coordinate in result
            w(i).d = 'W';                                                       %       Store wall direction in result
            w(i).p = p;                                                         %       Store probability in result
            w(i).color = obj.sco.walls(y,x,1).FaceVertexCData(1,:);             %       Store wall color in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
          p = obj.isWallAt(x,y,'S');                                            %     Get south wall prob. at location
          if p                                                                  %     If there is a west wall >>
            w(i).x = x;                                                         %       Store x coordinate in result
            w(i).y = y;                                                         %       Store y coordinate in result
            w(i).d = 'S';                                                       %       Store wall direction in result
            w(i).p = p;                                                         %       Store probability in result
            w(i).color = obj.sco.walls(y,x,2).FaceVertexCData(1,:);             %       Store wall color in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      w = w(1:i-1);                                                             % Shrink result
    end

    function removeAllWalls(obj)
      % Removes all walls.
      % 
      %   obj.removeAllWalls()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      for x=1:xdim+1                                                            % Loop over x coordinates >>
        for y=1:ydim+1                                                          %     Loop over y coordinates >>
          obj.removeWall(x,y,'W');                                              %     Remove west wall
          obj.removeWall(x,y,'S');                                              %     Remove east wall
        end                                                                     %   <<
      end                                                                       % <<
    end
    
    % - Mice

    function a=addMouse(obj,x,y,p,h)
      % Places a mouse at a location in the maze. If there is already a mouse 
      % at the specified location, the method will replace it. If the specified 
      % location is greater than the current maze limits in  either direction, 
      % the method will enlarge the maze.
      % 
      %   a=obj.addMouse(x,y)
      %   a=obj.addMouse(x,y,p)
      %   a=obj.addMouse(x,y,p,h)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   p   - Probability of mouse, 0 to remove, default is 1.
      %   h   - Heading angle of mouse in degrees, default is 0 (facing north).
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      if x<=0; error('x coordinate must be positive.'); end                     % Check x coordinate
      if y<=0; error('y coordinate must be positive.'); end                     % Check y coordinate
      if nargin<4; p = 1; end                                                   % Default probabilty: 1
      if nargin<5; h = 0; end                                                   % Default heading angle: 0°
      if obj.isMouseAt(x,y)                                                     % There is already a mouse at (x,y) >>
        delete(obj.sco.mice(y,x));                                              %   Delete from scene
        obj.sco.mice(y,x)=gobjects(1,1);                                        %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove tile
      if p>1; p=1; end                                                          % Limit probability to 1
      dx = 0; dy = 0; dz = 0;                                                   % Position offset of mouse
      if obj.isCheeseAt(x,y)                                                    % Cheese at new mouse position >>
        c = MouseMazePlot.getCheeseModel();                                     %   Get copy of cheese model
        dz = max(c.v); dz=dz(2);                                                %   Elevate mouse
      end                                                                       % <<
      if obj.isCheeseAt(x,y)||obj.isWaterAt(x,y)                                % Water or cheese at new mouse pos. >>
        dx = 0.08; dy = 0.05;                                                   %   Offset mouse on x-y plane
        h = obj.getMouseCamAngle();                                             %   Have mouse face camera
      end                                                                       % <<
      m = MouseMazePlot.getMouseModel();                                        % Get copy of mouse model
      m.v = m.v * MouseMazePlot.rotmat(h);                                      % Rotate
      m.v = m.v + [y+dy dz x+dx];                                               % Translate to maze coordinates
      a = patch(obj.ax,'vertices',m.v,'faces',m.f.v,'FaceVertexCData',m.fvcd);  % Create scene object
      f = fields(m.p); for i=1:length(f); a.set(f{i},m.p.(f{i})); end           % Copy patch properties
      a.FaceAlpha = p;                                                          % Probability -> face alpha
      obj.sco.mice(y,x) = a;                                                    % Store scene object
      obj.updateLabels();                                                       % Update row and column labels
    end

    function a=moveMouse(obj,x,y,h)
      % Moves mouse to new location in maze and removes all other mice.
      % 
      %   obj.moveMouse(x,y)
      %   obj.moveMouse(x,y,h)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   h   - Heading angle of mouse in degrees, default is 0 (facing north).
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      
      obj.removeAllMice();
      if nargin<4; h = 0; end                                                   % Default heading angle: 0°
      a = obj.addMouse(x,y,1,h);
    end

    function removeMouse(obj,x,y)
      % Removes a mouse. If there is no mouse at the specified location, the 
      % method will do nothing.
      %
      %   obj.removeMouse(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addMouse(x,y,0);                                                      % Invoke addMouse with p=0
    end
    
    function p=isMouseAt(obj,x,y)
      % Tests if there is a mouse at the specified location. If there is no 
      % mouse or if the specified location is outside the maze boundaries, the 
      % method will return 0.
      %
      %   p = obj.isMouseAt(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % returns:
      %   p   - The probability of a mouse at the specified location.
      p = 0;                                                                    % Default return value: 0
      if x<=0 || y<=0 || x>size(obj.sco.mice,2) || y>size(obj.sco.mice,1)       % Coords. outside scene object array >>
        return;                                                                 %   No cheese at (x,y)
      end                                                                       % <<
      if isa(obj.sco.mice(y,x),'matlab.graphics.primitive.Patch')               % Mouse patch at (x,y) >>
        p = obj.sco.mice(y,x).FaceAlpha;                                        %   Probability is alpha value
      end                                                                       % <<
    end

    function m=getAllMice(obj)
      % Gets all mice.
      %
      %   m=obj.getAllMice()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   m   - A struct vector containing the locations m(i).x, m(i).y and the
      %         probabilities m(i).p of all mice in the maze.
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = xdim*ydim;                                                            % Maximal number of mice
      m = struct('x', cell(1,n), 'y', cell(1,n), 'p', cell(1,n));               % Pre-allocate result
      i = 1;                                                                    % Pointer into result
      for x=1:xdim                                                              % Loop over x coordinates >>
        for y=1:ydim                                                            %   Loop over y coordinates >>
          p = obj.isMouseAt(x,y);                                               %     Get mouse prob. at location
          if p                                                                  %     If there is a mouse >>
            m(i).x = x;                                                         %       Store x coordinate in result
            m(i).y = y;                                                         %       Store y coordinate in result
            m(i).p = p;                                                         %       Store probability in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      m = m(1:i-1);                                                             % Shrink result
    end

    function removeAllMice(obj)
      % Removes all mice.
      % 
      %   obj.removeAllMice()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      for i=1:size(obj.sco.mice,2)
        for j=1:size(obj.sco.mice,1)
          obj.addMouse(i,j,0);
        end
      end
    end
    
    % - Cats

    function a=addCat(obj,x,y,p,h)
      % Places a cat at a location in the maze. If there is already a cat at 
      % the specified location, the method will replace it. If the specified 
      % location is greater than the current maze limits in  either direction, 
      % the method will enlarge the maze.
      % 
      %   a=obj.addCat(x,y)
      %   a=obj.addCat(x,y,p)
      %   a=obj.addCat(x,y,p,h)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   p   - Probability of cat, 0 to remove, default is 1.
      %   h   - Heading angle of cat in degrees, default is 0 (facing north).
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      if x<=0; error('x coordinate must be positive.'); end                     % Check x coordinate
      if y<=0; error('y coordinate must be positive.'); end                     % Check y coordinate
      if nargin<4; p = 1; end                                                   % Default probabilty: 1
      if nargin<5; h = 0; end                                                   % Default heading angle: 0°
      if obj.isCatAt(x,y)                                                       % There is already a cat at (x,y) >>
        delete(obj.sco.cats(y,x));                                              %   Delete from scene
        obj.sco.cats(y,x)=gobjects(1,1);                                        %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove tile
      if p>1; p=1; end                                                          % Limit probability to 1
      m = MouseMazePlot.getCatModel();                                          % Get copy of cat model
      m.v = m.v * MouseMazePlot.rotmat(h);                                      % Rotate
      m.v = m.v + [y 0 x];                                                      % Translate to maze coordinates
      a = patch(obj.ax,'vertices',m.v,'faces',m.f.v,'FaceVertexCData',m.fvcd);  % Create scene object
      f = fields(m.p); for i=1:length(f); a.set(f{i},m.p.(f{i})); end           % Copy patch properties
      a.FaceAlpha = p;                                                          % Probability -> face alpha
      obj.sco.cats(y,x) = a;                                                    % Store scene object
      obj.updateLabels();                                                       % Update row and column labels
    end

    function a=moveCat(obj,x,y,h)
      % Moves the cat to new location in maze and removes all other cats.
      % 
      %   obj.moveCat(x,y)
      %   obj.moveCat(x,y,h)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   h   - Heading angle of cat in degrees, default is 0 (facing north).
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      for i=1:size(obj.sco.cats,2)
        for j=1:size(obj.sco.cats,1)
          obj.addCat(i,j,0);
        end
      end
      if nargin<4; h = 0; end                                                   % Default heading angle: 0°
      a = obj.addCat(x,y,1,h);
    end

    function removeCat(obj,x,y)
      % Removes a cat. If there is no cat at the specified location, the method 
      % will do nothing.
      %
      %   obj.removeMouse(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addCat(x,y,0);                                                        % Invoke addCat with p=0
    end
    
    function p=isCatAt(obj,x,y)
      % Tests if there is a cat at the specified location. If there is no cat 
      % or if the specified location is outside the maze boundaries, the method 
      % will return 0.
      %
      %   p = obj.isCatAt(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % returns:
      %   p   - The probability of a cat at the specified location.
      p = 0;                                                                    % Default return value: 0
      if x<=0 || y<=0 || x>size(obj.sco.cats,2) || y>size(obj.sco.cats,1)       % Coords. outside scene object array >>
        return;                                                                 %   No cheese at (x,y)
      end                                                                       % <<
      if isa(obj.sco.cats(y,x),'matlab.graphics.primitive.Patch')               % Cat patch at (x,y) >>
        p = obj.sco.cats(y,x).FaceAlpha;                                        %   Probability is alpha value
      end                                                                       % <<
    end

    function c=getAllCats(obj)
      % Gets all cats.
      %
      %   c=obj.getAllCats()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   c   - A struct vector containing the locations c(i).x, c(i).y and the
      %         probabilities c(i).p of all mice in the maze.
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = xdim*ydim;                                                            % Maximal number of cats
      c = struct('x', cell(1,n), 'y', cell(1,n), 'p', cell(1,n));               % Pre-allocate result
      i = 1;                                                                    % Pointer into result
      for x=1:xdim                                                              % Loop over x coordinates >>
        for y=1:ydim                                                            %   Loop over y coordinates >>
          p = obj.isCatAt(x,y);                                                 %     Get cat probability at location
          if p                                                                  %     If there is a cat >>
            c(i).x = x;                                                         %       Store x coordinate in result
            c(i).y = y;                                                         %       Store y coordinate in result
            c(i).p = p;                                                         %       Store probability in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      c = c(1:i-1);                                                             % Shrink result
    end
    
    function removeAllCats(obj)
      % Removes all cats.
      % 
      %   obj.removeAllCats()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      for i=1:size(obj.sco.cats,2)
        for j=1:size(obj.sco.cats,1)
          obj.addCat(i,j,0);
        end
      end
    end

    % - Cheese
    
    function a=addCheese(obj,x,y,p)
      % Places cheese at a location in the maze.  If there is already
      % cheese at the specified location, the method will replace it. If the 
      % specified location is greater than the current maze limits in  either 
      % direction, the method will enlarge the maze.
      % 
      %   a=obj.addCheese(x,y)
      %   a=obj.addCheese(x,y,p)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   p   - Probability of cheese, 0 to remove, default is 1.
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      if x<=0; error('x coordinate must be positive.'); end                     % Check x coordinate
      if y<=0; error('y coordinate must be positive.'); end                     % Check y coordinate
      if nargin<4; p = 1; end                                                   % Default probabilty: 1
      if obj.isCheeseAt(x,y)                                                    % There is already cheese >>
        delete(obj.sco.cheese(y,x));                                            %   Delete from scene
        obj.sco.cheese(y,x)=gobjects(1,1);                                      %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove tile
      if p>1; p=1; end                                                          % Limit probability to 1
      c =  MouseMazePlot.getCheeseModel();                                      % Get copy of cheese model
      c.v = c.v+[y 0 x];                                                        % Translate to maze coordinates
      a = patch(obj.ax,'vertices',c.v,'faces',c.f.v,'FaceVertexCData',c.fvcd);  % Create new cheese scene object
      f = fields(c.p); for i=1:length(f); a.set(f{i},c.p.(f{i})); end           % Copy patch properties
      a.FaceAlpha = p;                                                          % Probability -> face alpha
      obj.sco.cheese(y,x) = a;                                                  % Store scene object
      obj.updateLabels();                                                       % Update row and column labels
    end

    function removeCheese(obj,x,y)
      % Removes cheese. If there is no cheese at the specified location, the 
      % method will do nothing.
      %
      %   obj.removeCheese(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addCheese(x,y,0);                                                     % Invoke addCheese with p=0
    end

    function p=isCheeseAt(obj,x,y)
      % Tests if there is cheese at the specified location. If there is no 
      % cheese or if the specified location is outside the maze boundaries, the 
      % method will return 0.
      %
      %   p = obj.isCheeseAt(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % returns:
      %   p   - The probability of cheese at the specified location.
      p = 0;                                                                    % Default return value: 0
      if x<=0 || y<=0 || x>size(obj.sco.cheese,2) || y>size(obj.sco.cheese,1)   % Coords. outside scene object array >>
        return;                                                                 %   No cheese at (x,y)
      end                                                                       % <<
      if isa(obj.sco.cheese(y,x),'matlab.graphics.primitive.Patch')             % Cheese patch at (x,y) >>
        p = obj.sco.cheese(y,x).FaceAlpha;                                      %   Probability is alpha value
      end                                                                       % <<
    end

    function c=getAllCheese(obj)
      % Gets all cheese bits.
      %
      %   c=obj.getCheese()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   c   - A struct vector containing the locations c(i).x, c(i).y and the
      %         probabilities c(i).p of all cheese bits in the maze.
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = xdim*ydim;                                                            % Maximal number of cheese bits
      c = struct('x', cell(1,n), 'y', cell(1,n), 'p', cell(1,n));               % Pre-allocate result
      i = 1;                                                                    % Pointer into result
      for x=1:xdim                                                              % Loop over x coordinates >>
        for y=1:ydim                                                            %   Loop over y coordinates >>
          p = obj.isCheeseAt(x,y);                                              %     Get cheese prob. at location
          if p                                                                  %     If there is a mouse >>
            c(i).x = x;                                                         %       Store x coordinate in result
            c(i).y = y;                                                         %       Store y coordinate in result
            c(i).p = p;                                                         %       Store probability in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      c = c(1:i-1);                                                             % Shrink result
    end

    function removeAllCheese(obj)
      % Removes all cheese.
      % 
      %   obj.removeAllCheese()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      for i=1:size(obj.sco.cheese,2)
        for j=1:size(obj.sco.cheese,1)
          obj.addCheese(i,j,0);
        end
      end
    end

    % - Water

    function a=addWater(obj,x,y,p)
      % Places water at a location in the maze.  If there is already water at 
      % the specified location, the method will replace it. If the specified 
      % location is greater than the current maze limits in either direction, 
      % the method will enlarge the maze.
      % 
      %   a=obj.addWater(x,y)
      %   a=obj.addWater(x,y,p)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %   p   - Probability of water, 0 to remove, default is 1.
      %
      % returns:
      %   a   - The newly created patch object;
      %
      % throws exception:
      % - If either x or y is non-positive.
      if x<=0; error('x coordinate must be positive.'); end                     % Check x coordinate
      if y<=0; error('y coordinate must be positive.'); end                     % Check y coordinate
      if nargin<4; p = 1; end                                                   % Default probabilty: 1
      if obj.isWaterAt(x,y)                                                     % There is already water >>
        delete(obj.sco.water(y,x));                                             %   Delete from scene
        obj.sco.water(y,x)=gobjects(1,1);                                       %   Remove from scene objects array
      end                                                                       % <<
      if p<=0; return; end                                                      % Just remove tile
      if p>1; p=1; end                                                          % Limit probability to 1
      c =  MouseMazePlot.getWaterModel();                                       % Get copy of water model
      c.v = c.v+[y 0 x];                                                        % Translate to maze coordinates
      a = patch(obj.ax,'vertices',c.v,'faces',c.f.v,'FaceVertexCData',c.fvcd);  % Create new cheese scene object
      f = fields(c.p); for i=1:length(f); a.set(f{i},c.p.(f{i})); end           % Copy patch properties
      a.FaceVertexAlphaData = p *a.FaceVertexAlphaData;                         % Probability -> face alpha data
      obj.sco.water(y,x) = a;                                                   % Store scene object
      obj.updateLabels();                                                       % Update row and column labels
    end

    function removeWater(obj,x,y)
      % Removes water. If there is no water at the specified location, the 
      % method will do nothing.
      %
      %   obj.removeWater(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.addWater(x,y,0);                                                      % Invoke addWater with p=0
    end

    function p=isWaterAt(obj,x,y)
      % Tests if there is water at the specified location. If there is no water 
      % or if the specified location is outside the maze boundaries, the method 
      % will return 0.
      %
      %   p = obj.isWaterAt(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % returns:
      %   p   - The probability of water at the specified location.
      p = 0;                                                                    % Default return value: 0
      if x<=0 || y<=0 || x>size(obj.sco.water,2) || y>size(obj.sco.water,1)     % Coords. outside scene object array >>
        return;                                                                 %   No water at (x,y)
      end                                                                       % <<
      if isa(obj.sco.water(y,x),'matlab.graphics.primitive.Patch')              % Water patch at (x,y) >>
        a = MouseMazePlot.getWaterModel();                                      %   Get water model
        ng3 = a.g{3}{2};                                                        %   Get first vertex of cup group
        p = obj.sco.water(y,x).FaceVertexAlphaData(ng3);                        %   Prob. is alpha value of that vertex
      end                                                                       % <<
    end

    function w=getAllWater(obj)
      % Gets all water cups.
      %
      %   c=obj.getWater()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %
      % returns:
      %   w   - A struct vector containing the locations w(i).x, w(i).y and the
      %         probabilities w(i).p of all water cups in the maze.
      [xdim,ydim] = obj.getDim();                                               % Get maze dimensions
      n = xdim*ydim;                                                            % Maximal number of water cups
      w = struct('x', cell(1,n), 'y', cell(1,n), 'p', cell(1,n));               % Pre-allocate result
      i = 1;                                                                    % Pointer into result
      for x=1:xdim                                                              % Loop over x coordinates >>
        for y=1:ydim                                                            %   Loop over y coordinates >>
          p = obj.isWaterAt(x,y);                                               %     Get water prob. at location
          if p                                                                  %     If there is a mouse >>
            w(i).x = x;                                                         %       Store x coordinate in result
            w(i).y = y;                                                         %       Store y coordinate in result
            w(i).p = p;                                                         %       Store probability in result
            i = i+1;                                                            %       Increment result pointer
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      w = w(1:i-1);                                                             % Shrink result
    end

    function removeAllWater(obj)
      % Removes all water.
      % 
      %   obj.removeAllWater()
      %
      % arguments:
      %   obj - The mouse-maze plot.
      
      for i=1:size(obj.sco.water,2)
        for j=1:size(obj.sco.water,1)
          obj.addWater(i,j,0);
        end
      end
    end    

  end

  %% == Graphics objects workers ==
  
  methods(Access=protected)
    
    function updateLabels(obj)
      % Updates the x and y place labels.
      [xdim,ydim]=obj.getDim();                                                 % Get maze dimensions
      if size(obj.sco.xlabels,1)~=xdim                                          % x dim. does not equal label count >>
        for x=1:size(obj.sco.xlabels,2); delete(obj.sco.xlabels(1,x)); end      %   Remove all x labels
        obj.sco.xlabels = gobjects(1,xdim);                                     %   Pre-alloc. x label scene object array
        for x=1:xdim                                                            %   Loop over x positions >>
          t = text(obj.ax,0.43,-0.07,x,sprintf('$|%d_X\\rangle$',x));           %     Create label
          t.Interpreter         = 'latex';                                      %     Interpret label text by LaTeX
          t.VerticalAlignment   = 'top';                                        %     Top-aligned
          t.HorizontalAlignment = 'center';                                     %     Centered
          obj.sco.xlabels(x)    = t;                                            %     Store scene object
        end                                                                     %   <<
      end                                                                       % <<
      if size(obj.sco.ylabels,1)~=ydim                                          % y dim. does not equal label count >>
        for y=1:size(obj.sco.ylabels,2); delete(obj.sco.ylabels(1,y)); end      %   Remove all y labels
        obj.sco.ylabels = gobjects(1,ydim);                                     %   Pre-alloc. x label scene object array
        for y=1:ydim                                                            %   Loop over y positions >>
          t = text(obj.ax,y,-0.07,0.43,sprintf('$|%d_Y\\rangle$',y));           %     Create label
          t.Interpreter         = 'latex';                                      %     Interpret label text by LaTeX
          t.VerticalAlignment   = 'middle';                                     %     Middle-aligned
          t.HorizontalAlignment = 'right';                                      %     Right-aligned
          obj.sco.ylabels(y)    = t;                                            %     Store scene object
        end                                                                     %   <<
      end                                                                       % <<
    end
    
    function alpha=getMouseCamAngle(obj)
      % Returns the rotation angle for the mouse facing the camera.
      %
      %   alpha=obj.getMouseCamAngle(obj)
      %
      % arguments:
      %   obj   - The mouse-maze plot.
      %
      % return:
      %   alpha - The mouse rotation angle in degrees.
      m = size(obj.sco.cheese)/2; % HACK: Actual mouse postion unknown.         % Get mouse position in maze coords.
      c = campos(obj.ax); c=[c(1,3) c(1,1)];                                    % Get camera position in maze coords.
      mc = c-m;                                                                 % Get mouse-to-camera vector
      mn = [0 1];                                                               % Get mouse-to-north vector
      alpha = acos(mn*mc'/(sqrt(mn*mn')*sqrt(mc*mc')))/pi*180;                  % Compute angle in degrees
      if mc(1,1)<0; alpha=360-alpha; end                                        % Correct if camera left of mouse
    end
    
  end

  methods(Static)%, Access=protected)

    function obj=getFloorTileModel()
      % Retrieves the floor tile object.
      persistent floorTileModel;                                                % Persistent floor tile model variable
      if isempty(floorTileModel)                                                % Not yet initialized >>
        tt = 0.03;                                                              %   Tile thickness
        floorTileModel.v = [                                                    %   Vertices >>
          -0.5+tt/2 0     -0.5+tt/2;                                            %      #1 (face, bottom-left)
           0.5-tt/2 0     -0.5+tt/2;                                            %      #2 (face, bottom-right)
           0.5-tt/2 0      0.5-tt/2;                                            %      #3 (face, top-right)
          -0.5+tt/2 0      0.5-tt/2;                                            %      #4 (face, top-left)
          -0.5      -tt/2 -0.5;                                                 %      #5 (middle, bottom-left)
           0.5      -tt/2 -0.5;                                                 %      #6 (middle, bottom-right)
           0.5      -tt/2  0.5;                                                 %      #7 (middle, top-right)
          -0.5      -tt/2  0.5;                                                 %      #8 (middle, top-left)
          -0.5      -tt   -0.5;                                                 %      #9 (back, bottom-left)
           0.5      -tt   -0.5;                                                 %     #10 (back, bottom-right)
           0.5      -tt    0.5;                                                 %     #11 (back, top-right)
          -0.5      -tt    0.5                                                  %     #12 (back, top-left)
        ];                                                                      %   <<
        floorTileModel.f.v = [                                                  %   Faces >>
           1  2  3  4;                                                          %     Face
           5  6  2  1;                                                          %     Bevel (bottom)
           6  7  3  2;                                                          %     Bevel (right)
           7  8  4  3;                                                          %     Bevel (top)
           8  5  1  4;                                                          %     Bevel (left)
           9 10  6  5;                                                          %     Side (bottom)
          10 11  7  6;                                                          %     Side (right)
          11 12  8  7;                                                          %     Side (top)
          12  9  5  8;                                                          %     Side (left)
           9 10 11 12                                                           %     Back
        ];                                                                      %   <<
        fclr = MouseMazePlot.clrFloor+0.75*(1-MouseMazePlot.clrFloor);          %   Make face color
        sclr = 0.3.*MouseMazePlot.clrFloor;                                     %   Make sides color
        bclr = fclr;                                                            %   Make back color
        floorTileModel.fvcd = repmat(fclr,4,1);                                 %   Set face color
        floorTileModel.fvcd = [floorTileModel.fvcd; repmat(sclr,4,1)];          %   Set bevel color
        floorTileModel.fvcd = [floorTileModel.fvcd; repmat(bclr,4,1)];          %   Set side and back color
        floorTileModel.p.FaceColor        = 'flat';                             %   Flat face colors
        floorTileModel.p.FaceLighting     = 'gouraud';                          %   Flat face shading
        floorTileModel.p.EdgeColor        = sclr;                               %   Show edges (uniform color)
        floorTileModel.p.EdgeLighting     = 'flat';                             %   Uniform edge shading
        floorTileModel.p.LineWidth        = 0.01;                               %   Edge hairlines
        floorTileModel.p.DiffuseStrength  = 0.9;                                %   Strength of diffuse light
        floorTileModel.p.AmbientStrength  = 0.8;                                %   Strength of ambient light
        floorTileModel.p.SpecularStrength = 0;                                  %   Dull
      end                                                                       % <<
      obj = floorTileModel;                                                     % Return copy of persistent variable
    end
    
    function obj=getWallModel()
      persistent wallModel;                                                     % Persistent wall model variable
      if isempty(wallModel)                                                     % Not yet initialized >>
        ft = MouseMazePlot.getFloorTileModel();                                 %   Get floor tile model
        tt = max(ft.v)-min(ft.v); tt = tt(2);                                   %   Get thickness of floor tile
        wt = 3*tt;                                                              %   Wall thickness: 2x tile thickness
        wh = 0.6;                                                               %   Wall height
        wallModel.v = [                                                         %   Vertices >>
          -wt/2  0   0.5;                                                       %      #1
          -wt/2  0  -0.5;                                                       %      #2
          -wt/2  wh -0.5;                                                       %      #3
          -wt/2  wh  0.5;                                                       %      #4
           wt/2  0   0.5;                                                       %      #5
           wt/2  0  -0.5;                                                       %      #6
           wt/2  wh -0.5;                                                       %      #7
           wt/2  wh  0.5;                                                       %      #8
          -wt/2  0   0.5;                                                       %      #9 (Repetition of 1-8)
          -wt/2  0  -0.5;                                                       %     #10
          -wt/2  wh -0.5;                                                       %     #11
          -wt/2  wh  0.5;                                                       %     #12
           wt/2  0   0.5;                                                       %     #13
           wt/2  0  -0.5;                                                       %     #14
           wt/2  wh -0.5;                                                       %     #15
           wt/2  wh  0.5                                                        %     #16
        ];                                                                      %   <<
        wallModel.f.v = [                                                       %   Faces >>
           1  2  3  4;                                                          %     Left face
           5  6  7  8;                                                          %     Right face
           9 10 14 13;                                                          %     Bottom side
          10 14 15 11;                                                          %     Front side
          11 15 16 12;                                                          %     Top side
          12 16 13  9;                                                          %     Rear side
        ];                                                                      %   <<
        fclr = MouseMazePlot.clrFloor+0.75*(1-MouseMazePlot.clrFloor);          %   Make face color
        sclr = 0.3.*MouseMazePlot.clrFloor;                                     %   Make sides color
        wallModel.fvcd = repmat(fclr,8,1);                                      %   Set faces color
        wallModel.fvcd = [wallModel.fvcd; repmat(sclr,8,1)];                    %   Set sides color
        wallModel.p.FaceColor        = 'flat';                                  %   Flat face colors
        wallModel.p.FaceLighting     = 'gouraud';                               %   Flat face shading
        wallModel.p.EdgeColor        = sclr;                                    %   Show edges (uniform color)
        wallModel.p.EdgeLighting     = 'flat';                                  %   Uniform edge shading
        wallModel.p.LineWidth        = 0.01;                                    %   Edge hairlines
        wallModel.p.DiffuseStrength  = 0.9;                                     %   Strength of diffuse light
        wallModel.p.AmbientStrength  = 0.8;                                     %   Strength of ambient light
        wallModel.p.SpecularStrength = 0;                                       %   Dull
      end                                                                       % <<
      obj = wallModel;                                                          % Return copy of persistent variable
    end
    
    function obj=getMouseModel() 
      % Retrieves the mouse 3D model.
      persistent mouseModel;                                                    % Persistent mouse model variable
      if isempty(mouseModel)                                                    % Not yet initialized >>
        mouseModel = MouseMazePlot.fromCache('mouse');                          %   Try to get mouse object from cache
        if isempty(mouseModel)                                                  %   Not in chace >>
          fprintf('Loading and caching mouse model ...');                       %     Screen log
          scale = 0.004;                                                        %     Vertex scale factor
          mouseModel = readObj('mouse.obj','mouse.png');                        %     Load mouse object
          mouseModel.v = mouseModel.v * [scale 0 0; 0 scale 0; 0 0 scale];      %     Scale mouse object
          mouseModel.v = mouseModel.v * MouseMazePlot.rotmat(90);               %     Rotate mouse object to face north
          mouseModel.v = mouseModel.v + [-0.25 0 0];                            %     Translate mouse object
          mouseModel.p.FaceColor        = 'interp';                             %     Interpolate face colors
          mouseModel.p.FaceLighting     = 'gouraud';                            %     Gouraud face shading
          mouseModel.p.EdgeColor        = 'none';                               %     Hide edges
          mouseModel.p.DiffuseStrength  = 0.9;                                  %     Strength of diffuse light
          mouseModel.p.AmbientStrength  = 0.3;                                  %     Strength of ambient light
          mouseModel.p.SpecularStrength = 0.2;                                  %     Weak shine
          MouseMazePlot.toCache('mouse',mouseModel);                            %     Cache mouse model
          fprintf(' done\n');                                                   %     Screen log
        end                                                                     %   <<
      end                                                                       % <<
      obj = mouseModel;                                                         % Return copy of persistent variable
    end

    function obj=getCatModel()
      % Retrieves the cat 3D model.
      persistent catModel;                                                      % Persistent cat model variable
      if isempty(catModel)                                                      % Not yet initialized >>
        catModel = MouseMazePlot.fromCache('cat');                              %   Try to get cat object from cache
        if isempty(catModel)                                                    %   Not in cache >>
          fprintf('Loading and caching cat model ...');                         %     Screen log
          scale = 0.0018;                                                       %     Vertex scale factor
          catModel = readObj('cat.obj','cat.png');                              %     Load cat object
          catModel.v = catModel.v * [0 0 scale; 0 scale 0; scale 0 0];          %     Scale cat object
          catModel.v = catModel.v + [0.15 0 0];                                 %     Translate cat object
          catModel.p.FaceColor        = 'interp';                               %     Interpolate face colors
          catModel.p.FaceLighting     = 'gouraud';                              %     Gouraud face shading
          catModel.p.EdgeColor        = 'none';                                 %     Hide edges
          catModel.p.DiffuseStrength  = 1;                                      %     Strength of diffuse light
          catModel.p.AmbientStrength  = 0.8;                                    %     Strength of ambient light
          catModel.p.SpecularStrength = 0.2;                                    %     Weak shine
          MouseMazePlot.toCache('cat',catModel);                                %     Cache object
          fprintf(' done\n');                                                   %     Screen log
        end                                                                     %   <<
      end                                                                       % <<
      obj = catModel;                                                           % Return copy of persistent variable
    end

    function obj=getCheeseModel()
      % Retrieves the cheese 3D model.
      persistent cheeseModel;                                                   % Persistent cheese model variable
      if isempty(cheeseModel)                                                   % Not yet initialized >>
        cheeseModel = MouseMazePlot.fromCache('cheese');                        %   Try to get cheese object from cache
        if isempty(cheeseModel)                                                 %   Not in cache >>
          fprintf('Loading and caching cheese model ...');                      %     Screen log
          scale = 0.1;                                                          %     Vertex scale factor
          cheeseModel = readObj('cheese.obj',[1 0.8 0]);                        %     Load cheese object
          cheeseModel.v = cheeseModel.v * [scale 0 0; 0 scale 0; 0 0 scale];    %     Scale cheese object
          mn = min(cheeseModel.v);                                              %     Get minimum of vertex coordinates
          cheeseModel.v = cheeseModel.v + [0 -mn(2) 0.1];                       %     Translate object
          cheeseModel.p.FaceColor       = 'interp';                             %     Interpolate face colors
          cheeseModel.p.FaceLighting    = 'gouraud';                            %     Gouraud face shading
          cheeseModel.p.EdgeColor       = 'none';                               %     Hide edges
          cheeseModel.p.AmbientStrength = 0.35;                                 %     Strength of ambient light
          MouseMazePlot.toCache('cheese',cheeseModel);                          %     Cache object
          fprintf(' done\n');                                                   %     Screen log
        end                                                                     %   <<
      end                                                                       % <<
      obj = cheeseModel;                                                        % Return copy of persistent variable
    end

    function obj=getWaterModel()
      % Retrieves the water glass 3D model.
      persistent waterModel;                                                    % Persistent waterglass model variable
      if isempty(waterModel)                                                    % Not yet initialized >>
        waterModel = MouseMazePlot.fromCache('water');                          %   Try to get watergl.object from cache
        if isempty(waterModel)                                                  %   Not in cache >>
          fprintf('Loading and caching waterglass model ...');                  %     Screen log
          scale = 0.00025;                                                      %     Vertex scale factor
          waterModel = readObj('water.obj',[0.7 0.9 1]);                        %     Load waterglass object
          waterModel.v = waterModel.v * [scale 0 0; 0 scale 0; 0 0 scale];      %     Scale waterglass object
          waterModel.v = waterModel.v + [-0.1 0 -0.1];                          %     Translate object
          waterModel.p.FaceColor        = 'interp';                             %     Interpolate face colors
          waterModel.p.FaceAlpha        = 'interp';                             %     Interpolate face transparencies
          waterModel.p.AlphaDataMapping = 'none';                               %     No alpha data mapping
          waterModel.p.FaceLighting     = 'gouraud';                            %     Gouraud face shading
          waterModel.p.EdgeColor        = 'none';                               %     Hide edges
          waterModel.p.AmbientStrength  = 0.4;                                  %     Strength of ambient light
          waterModel.p.DiffuseStrength  = 0.7;                                  %     Strength of diffuse light
          waterModel.p.SpecularStrength = 0.9;                                  %     Specular strength
          waterModel.p.SpecularExponent = 10;                                   %     Specular exponent
          waterModel.p.SpecularColorReflectance = 1;                            %     Specular color reflectance
          nv  = size(waterModel.v,1);                                           %     Get number of vertices
          ng3 = waterModel.g{3}{2};                                             %     Get first vertex of cup
          waterModel.p.FaceVertexAlphaData = ...                                %     Set alphas for water and cup
            [ 0.3*ones(ng3-1,1); ones(nv-ng3+1,1) ];                            %     ...     
          waterModel.fvcd(ng3:nv,:) = repmat([1 1 1],nv-ng3+1,1);               %     Set color of cup -> white
          MouseMazePlot.toCache('water',waterModel);                            %     Cache object
          fprintf(' done\n');                                                   %     Screen log
        end                                                                     %   <<
      end                                                                       % <<
      obj = waterModel;                                                         % Return copy of persistent variable
    end

    function rm=rotmat(alpha)
      % Computes a mouse rotation matrix.
      rm = [
        cos(alpha/180*pi)  0 sin(alpha/180*pi);
        0                  1 0;
        -sin(alpha/180*pi) 0 cos(alpha/180*pi)];
    end
    
    function [x,y,d]=convertWallCoords(x,y,d)
      assert(x>0&&y>0,'Maze coordinates must be positive');                     % Check x and y coordinates
      if     strcmpi(d,'W'); d = 1;                                             % West side
      elseif strcmpi(d,'S'); d = 2;                                             % South side
      elseif strcmpi(d,'E'); d = 1; x = x+1;                                    % West side -> west side of next place
      elseif strcmpi(d,'N'); d = 2; y = y+1;                                    % North side-> south side of next place
      else;  error('Argument d must be one of ''N'', ''E'', ''S'', or ''W''.'); % All other value of d invalid
      end                                                                       %
    end
    
  end
  
  %% == Patch cache ==

  methods(Static, Access=protected)
    
    function toCache(name,obj)
      assert(~isempty(obj));
      fname = [MouseMazePlot.cacheDir '/' name '.mat'];
      if exist(MouseMazePlot.cacheDir,'dir')~=7
        mkdir(MouseMazePlot.cacheDir);
      end
      save(fname,'obj');
    end
    
    function obj=fromCache(name)
      fname = [MouseMazePlot.cacheDir '/' name '.mat'];
      if exist(fname,'file')==2
        load(fname,'obj');
        return;
      else
        obj = [];
      end
    end
    
  end

  %% == Extended static API ==

  methods(Static)

    function clearCache()
      if exist(MouseMazePlot.cacheDir,'dir')==7
        rmdir(MouseMazePlot.cacheDir,'s');
      end
    end

    function a=preview(o)
      % Preview 3D object.
      figure;
      axis equal;
      xlabel('x'); ylabel('y'); zlabel('z');
      a = patch('vertices',o.v,'faces',o.f.v,'FaceVertexCData',o.fvcd);
      f = fields(o.p); for i=1:length(f); a.set(f{i},o.p.(f{i})); end
      camproj('perspective');
      camlight;
      %shading interp
      %lighting gouraud;
      cameratoolbar;
      cameratoolbar('SetCoordSys','y');
      
      camdolly(0,8,0,'fixtarget');
      camorbit(65,0,'data',[0 1 0]);
      % Make MP4 video -->
      % v = VideoWriter('mouse.mp4','MPEG-4');
      % open(v);
      % for i=1:360
      %   camorbit(1,0,'data',[0 1 0])
      %   drawnow
      %   frame = getframe(gcf);
      %   writeVideo(v,frame);
      % end
      % close(v);
      % <--
    end

    function r=DEBUG()
      % --Debugging use only--
      
      % Display cat texture grid
      xx = 512;
      yy = 512;
      r = MouseMazePlot.getCatModel();
      C = ones(xx+1,yy+1,3);
      for i=1:size(r.vt,1)
        x = max(1,floor(r.vt(i,1)*xx)+1);
        y = max(1,floor(r.vt(i,2)*yy)+1);
        fprintf('x=%d, y=%d\n',x,y);
        C(x,y,:)=[0 0 0];
      end
      image(C);
      imwrite(C,'cat_texturegrid.png');
    end

  end

end

% EOF