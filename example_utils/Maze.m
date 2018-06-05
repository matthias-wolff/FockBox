classdef Maze < handle
  % A maze world in which agents (mice or cats) live.
  %
  % See also: Agent
  
  %% == Properties ==
  properties(SetAccess=protected)
    
    % MouseMazePlot visualizing the world.
    mmp;

  end
    
  properties(Access=protected)
    
    % Map of agents living in the maze world.
    agents;
    
  end

  %% == Public API ==
  methods
 
    function obj=Maze()
      % Creates a maze world.
      
      obj.agents = containers.Map;
      figure('Name','World','NumberTitle','off');
      obj.mmp = MouseMazePlot([0 0]);
    end

    % TODO: Write documentation comments
    function o = addAgent(obj,agent,x,y)
      % TODO: ...
      if nargin < 4                                                             % No coordinates given >>
        x = 1;                                                                  %   Use default
        y = 1;                                                                  %   ...
      end                                                                       % <<
      agent.world = obj;                                                        % Set agent's world
      props.agent = agent;                                                      % Create agent's properties
      props.x = x;                                                              % ...
      props.y = y;                                                              % ...
      obj.agents(agent.id) = props;                                             % Store agent's properties
      o = obj.teleport(agent,x,y);                                              % Move agent
    end
    
    function addPlace(obj,x,y)
      % Adds a place to the maze world.
      % 
      %   obj.addPlace(x,y)
      %
      % arguments:
      %   obj - The mouse-maze plot.
      %   x   - x coordinate of maze location.
      %   y   - y coordinate of maze location.
      %
      % throws exception:
      % - If either x or y is non-positive.

      obj.mmp.addFloorTile(x,y);
    end

    function addWall(obj,x,y,d)
      % Adds a wall to the maze world.
      % 
      %   obj.addWall(x,y,d)
      %
      % arguments:
      %   obj   - The mouse-maze plot.
      %   x     - x coordinate of maze location.
      %   y     - y coordinate of maze location.
      %   d     - The wall position, one of 'N', 'E', 'S', or 'W'.
      %
      % throws exception:
      % - If either x or y is non-positive.
      obj.mmp.addWall(x,y,d,0.3);
    end

    % TODO: Write documentation comments
    function addCheese(obj,x,y)
      % TODO: ...
      obj.mmp.addCheese(x,y);
    end
    
    % TODO: Write documentation comments
    function addWater(obj,x,y)
      % TODO: ...
      obj.mmp.addWater(x,y);
    end

  end

  %% == API available to agents only ==
  methods(Access={?Agent})
    
    function o=action(obj,agent,a)
      % Performs an action in the maze.
      %
      %    o = obj.action(agent)
      %    o = obj.action(agent,a)
      %
      % arguments:
      %    obj   - The maze world.
      %    agent - The agent to perform the action with.
      %    a     - The action, a tensor product of a ket vector for the
      %            movement in x-direction (|+X>, |0X|, or |-X>) and a
      %            ket vector for the movement in y-direction (|+Y>, |0Y|,
      %            or |-Y>). The default is |0X0Y>.
      %
      % returns:
      %    o     - Observation after performing the action (a ket).
      %
      % example:
      %
      %    obj.action(agent,fockobj.bket('+X§0Y')); % Go east
      
      % Initialize                                                              % -------------------------------------
      obj.checkAgent(agent);                                                    % Check agent
      x  = obj.agents(agent.id).x;                                              % Get agent's current x-coordinate
      y  = obj.agents(agent.id).y;                                              % Get agent's current y-coordinate
      dx = 0;                                                                   % Movement delta on x axis
      dy = 0;                                                                   % Movement delta on y axis
      h  = 200;                                                                 % Agent's heading angle

      % Decode action                                                           % -------------------------------------
      if nargin>2                                                               % Action ket passed >>
        Xinc = fockobj.bket('+X');                                              %   Increment x basis vector
        Xnop = fockobj.bket('0X');                                              %   No op. on x basis vector
        Xdec = fockobj.bket('-X');                                              %   Decrement x basis vector
        Yinc = fockobj.bket('+Y');                                              %   Increment y basis vector
        Ynop = fockobj.bket('0Y');                                              %   No op. on y basis vector
        Ydec = fockobj.bket('-Y');                                              %   Decrement y basis vector
        Px   = Xinc*[Xinc (Yinc+Ynop+Ydec)]' ...                                %   x movement projector
             + Xnop*[Xnop (Yinc+Ynop+Ydec)]' ...                                %   ...
             + Xdec*[Xdec (Yinc+Ynop+Ydec)]';                                   %   ...
        Py   = Yinc*[(Xinc+Xnop+Xdec) Yinc]' ...                                %   y movement projector
             + Ynop*[(Xinc+Xnop+Xdec) Ynop]' ...                                %   ...
             + Ydec*[(Xinc+Xnop+Xdec) Ydec]';                                   %   ...
        if a'*Px'*Xinc~=0; dx=1; elseif a'*Px'*Xdec~=0; dx=-1; end              %   Get delta x of movement
        if a'*Py'*Yinc~=0; dy=1; elseif a'*Py'*Ydec~=0; dy=-1; end              %   Get delta y of movement
      end                                                                       % <<

      % Compute new coordinates of agent                                        % -------------------------------------
      if dx>0 && obj.mmp.isWallAt(x,y,'E'); dx = 0; end                         % Check for eastern wall
      if dx<0 && obj.mmp.isWallAt(x,y,'W'); dx = 0; end                         % Check for western wall
      if dy>0 && obj.mmp.isWallAt(x,y,'N'); dy = 0; end                         % Check for northern wall
      if dy<0 && obj.mmp.isWallAt(x,y,'S'); dy = 0; end                         % Check for southern wall
      xn = x + dx;                                                              % Agent's new x position
      yn = y + dy;                                                              % Agent's new y position
      if xn<1; xn = obj.mmp.getDim(1); end                                      % Wrap left
      if xn>obj.mmp.getDim(1); xn = 1; end                                      % Wrap right
      if yn<1; xn = obj.mmp.getDim(2); end                                      % Wrap bottom
      if yn>obj.mmp.getDim(2); yn = 1; end                                      % Wrap top
      if dx~=0 || dy~=0; h = angle(dx*1i+dy)/pi*180; end                        % Compute agent's heading angle

      % Perform action                                                          % -------------------------------------
      o = obj.teleport(agent,xn,yn,h);                                          % Move agent

    end

    % TODO: Write documentation comments
    function o=teleport(obj,agent,x,y,h)
      % TODO: Implement!

      % Initialize                                                              % -------------------------------------
      if nargin < 5; h = 200; end                                               % Default heading angle of agent
      obj.checkAgent(agent);                                                    % Check agent
      xn = x;                                                                   % New x-coordinate of agent
      yn = y;                                                                   % New y-coordinate of agent
      x  = obj.agents(agent.id).x;                                              % Get agent's current x-coordinate
      y  = obj.agents(agent.id).y;                                              % Get agent's current y-coordinate

      % Move agent                                                              % -------------------------------------
      props = obj.agents(agent.id);                                             % Get agent properties
      props.x = xn;                                                             % Set agent's new x position
      props.y = yn;                                                             % Set agent's new y position
      obj.agents(agent.id) = props;                                             % Set agent properties
      if agent.isMouse()                                                        % Agent is a mouse >>
        obj.mmp.removeMouse(x,y);                                               %   Move mouse agent in plot
        obj.mmp.addMouse(xn,yn,1,h);                                            %   ...
      else                                                                      % << TODO >>
        % TODO: Other agent types!
      end                                                                       % <<
      drawnow;                                                                  % Update plot

      % Make observation                                                        % -------------------------------------
      x = fockobj.bket(sprintf('%dX',xn));                                      % Make x-coordinate ket
      y = fockobj.bket(sprintf('%dY',yn));                                      % Make y-coordinate ket
      c = fockobj.bket('0C');                                                   % Make cheese ket
      if obj.mmp.isCheeseAt(xn,yn); c = fockobj.bket('1C'); end                 % ...
      g = fockobj.bket('0G');                                                   % Make water ket
      if obj.mmp.isWaterAt(xn,yn); g = fockobj.bket('1G'); end                  % ...
      o = [x y c g];                                                            % Make observaetion ket
    end
    
  end

  %% == Auxiliary methods ==
  methods(Access=protected)

    function checkAgent(obj,agent)
      % Checks an agent object.
      
      if ~isa(agent,'Agent')
        error('Argument is not an Agent object');
      end
      if ~obj.agents.isKey(agent.id)
        error('Agent is not in maze. Invoke addAgent first!');
      end
      agent.world = obj;
    end
  
  end
end

% EOF