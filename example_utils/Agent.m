classdef Agent < handle
  % An agent.
  
  properties(SetAccess=protected)
    
    % Unique identifier of agent.
    id;
    
    % MouseMazePlot visualizing the agent's inner stage.
    mmp;

    % The world state.
    w;
    
    % The agent state.
    phi;
    
    % The action repertoire.
    act;

    % Operator cache.
    Oc;
    
  end
  
  properties(SetAccess={?Maze})
    
    % The Maze world agent is living in.
    world;
    
  end
  
  % == Instance getters and information ==
  methods(Static)
    
    function obj=createMouse()
      % Creates a mouse agent.
      %
      %    obj = Agent.createMouse();
      %
      % returns:
      %    obj - The newly created agent.
      
      obj = Agent('mouse');
    end
    
  end
  
  % == Setters and getters ==
  methods

    function b=isMouse(obj)
      % Determines if an agent is a mouse.
      %
      %   b = obj.isMouse()
      %
      % returns:
      %   True if the agent is a mouse, false otherwise.
      
      b = startsWith(obj.id,'mouse');
    end
    
  end
  
  % == Perception-action cycle
  methods

    function addAction(obj,aid,a)
      % Adds an action to the agent's repertoire.
      %
      %   obj.addAction(aid,a)
      %
      % arguments:
      %   obj - The agent.
      %   aid - The action identifier.
      %   a   - The action semantics as understood by the world (a ket).
      
      action.a = a;                                                             % Create world action
      action.O = 0;                                                             % Create trial action
      obj.act(aid) = action;                                                    % Add action
    end

    function o=teleport(obj,x,y)
      % Teleports the agent to a location in the world.
      %
      %   o = obj.action(aid,x,y)
      %
      % arguments:
      %   obj - The agent.
      %   x   - x-coordinate of location to teleport to.
      %   y   - y-coordinate of location to teleport to.
      %
      % returns:
      %   o   - Observation of world in result of action (a ket).

      fprintf('___________________________________________________________\n'); % Console log
      fprintf('TELEPORT agent "%s" to (%d,%d):\n',obj.id,x,y);                  % Console log
      o = obj.perceive(obj.world.teleport(obj,x,y));                            % Teleport and preceive result
    end

    function o=action(obj,aid)
      % Performs an action in the world.
      %
      %   o = obj.action(aid)
      %
      % arguments:
      %   obj - The agent.
      %   aid - Identifier of action.
      %
      % returns:
      %   o   - Observation of world in result of action (a ket).

      a = obj.act(aid).a;                                                       % Get action semantics
      fprintf('___________________________________________________________\n'); % Console log
      fprintf('ACTION "%s" of agent "%s":\n',aid,obj.id);                       % Console log
      disp(a);                                                                  % Console log
      o = obj.perceive(obj.world.action(obj,a),aid);                            % Execute action and perceive result
    end

    function o=perceive(obj,o,aid)
      % Perceive an observation from the world.
      %
      %    o = obj.perceive(o,aid)
      %    o = obj.perceive(o)
      %
      % arguments:
      %    o   - The observation semantics as output by the world (a ket).
      %    aid - Identifier of action that caused the perception (optional).
      %
      % returns:
      %    o

      % Initialize                                                              % -------------------------------------
      pw = 0;                                                                   % Wall-perceived flag
      h = 200;                                                                  % Default agent heading

      % Update trial action operator                                            % -------------------------------------
      if nargin>=3                                                              % Action caused observation >>
        % TODO: Update trial action operator
        [x0,y0] = obj.locate();                                                 %   Locate agent before observation
        if o==obj.phi; pw = 0.3; end                                            %   Action ineffective -> wall
        if     aid=='N'; h = 0;                                                 %   Agent heading for north action
        elseif aid=='E'; h = 90;                                                %   Agent heading for east action
        elseif aid=='S'; h = 180;                                               %   Agent heading for south action
        elseif aid=='W'; h = -90;                                               %   Agent heading for west action
        end                                                                     %
      end                                                                       % <<

      % Update agent and world states                                           % -------------------------------------
      obj.phi = o;                                                              % Update agent state
      obj.w = obj.w + (1-obj.w'*o)*o;                                           % Update world state

      % Update plot                                                             % -------------------------------------
      if nargin>=3
        obj.renderInnerStage(aid);
      else
        obj.renderInnerStage();
      end
      
      % Perceive cheese and water                                               % -------------------------------------
      [x,y] = obj.locate();                                                     % Locate agent after observation
      pc = obj.phi' * obj.getLocationProjector(x,y,1,nan) * obj.phi;            % Probability of cheese at (x,y)
      pg = obj.phi' * obj.getLocationProjector(x,y,nan,1) * obj.phi;            % Probability of water at (x,y)
      
      % Update plot                                                             % -------------------------------------
      obj.mmp.addFloorTile(x,y);                                                % Add place to plot
      obj.mmp.addCheese(x,y,pc);                                                % Add/remove cheese
      obj.mmp.addWater(x,y,pg);                                                 % Add/remove water
      if obj.isMouse                                                            % Agent is a mouse >>
        obj.mmp.moveMouse(x,y,h);                                               %   Move mouse in plot
      else                                                                      % << TODO >>
        % TODO: other agent types
      end                                                                       % <<
      if nargin>=3; obj.mmp.addWall(x0,y0,aid,pw); end                          % Add/remove wall
      drawnow;                                                                  % Update plots
      fprintf('\n> Agent state:\n'); disp(obj.phi);                             % Console log
      fprintf('\n> World state:\n'); disp(obj.w);                               % Console log
    end

    function rexplore(obj,n)
      % Performs a random exploration.
      %
      %   obj.rexplore(n)
      %
      % arguments:
      %   obj - The agent.
      %   n   - The number of exploration actions.
      
      aids = obj.act.keys;                                                      % Get action ids
      for i=1:n                                                                 % Random exploration iterations >>
        aid = aids{randi(length(aids))};                                        %   Randomly select action id
        obj.action(aid);                                                        %   Perform action
      end                                                                       % <<
      
    end

    function sexplore(obj,rx,ry)
      % Performs a systematic exploration.
      %
      %   obj.sexplore(rx,ry)
      %
      % arguments:
      %   obj - The agent.
      %   rx  - Range of x-coordinates to explore, e.g. 1:10.
      %   ry  - Range of y-coordinates to explore.
      
      aids = obj.act.keys;                                                      % Get action ids
      for y=ry                                                                  % Loop over x-coordinates >>
        for x=rx                                                                %   Loop over y-coordinates >>
          for i=1:length(aids)                                                  %     Loop over actions >>
            obj.teleport(x,y);                                                  %       Teleport agent to coordinates
            obj.action(aids{i});                                                %       Perform action
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      
    end

  end
  
  % == Inner stage ==
  methods

    function p=isat(obj,x,y)
      % Determines the probability that the agent is on specified location.
      %
      %   p = obj.isat(x,y)
      %
      % arguments:
      %   obj - The agent.
      %   x   - The x-coordinate of queried location.
      %   y   - The y-coordinate of queried location.
      %
      % returns:
      %   The probability that the agent is on specified location.
      
      P = obj.getLocationProjector(x,y);                                        % Get location projector
      p = obj.phi' * P *obj.phi;                                                % Measure agent state
    end

    function [x,y,p]=locate(obj)
      % Determines the most probable location of the agent.
      %
      %   [x,y,p] = obj.locate()
      %
      % arguments:
      %   obj - The agent.
      %
      % returns:
      %   x   - The x-coordinate of the most probable location.
      %   y   - The y-coordinate of the most probable location.
      %   p   - The probability.
      
      xx = obj.mmp.getDim(1) + 1; yy = obj.mmp.getDim(2) + 1;                   % Get maze size
      x = 0; y = 0; p = -1;                                                     % Default return values
      for xi=1:xx                                                               % Loop over x-coordinates >>
        for yi=1:yy                                                             %   Loop over y-coordinates >>
          pi = obj.isat(xi,yi);                                                 %     Prob. of agent being @(xi,yi)
          if pi>p                                                               %     Greatest prob. so far >>
            x = xi; y = yi; p = pi;                                             %       Update return values
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
    end

    function P=getLocationProjector(obj,x,y,c,g)
      % Returns a location projector
      %
      %   P = obj.getLocationProjector(x,y,c,g)
      %
      % arguments:
      %   obj - The agent.
      %   x   - The x-coordinate of the location.
      %   y   - The y-coordinate of the location.
      %   c   - Cheese flag: 0 - no cheese, 1 - cheese, otherwise - dont' care
      %   g   - Water flag: 0 - no water, 1 - water, otherwise - dont' care
      %
      % returns:
      %   P   - The location projector.
      
      if nargin<4; c=-1; end; if c~=0&&c~=1; c=-1; end                          % Default cheese flag, normalize
      if nargin<5; g=-1; end; if g~=0&&g~=1; g=-1; end                          % Default water flag, normalize
      key = sprintf('PL(%d,%d,%d,%d)',x,y,c,g);                                 % Make operator cache key
      if obj.Oc.isKey(key)                                                      % Location projector in cache >>
        P = obj.Oc(key);                                                        %   Retrieve location projector
      else                                                                      % << Location projector not in cache >>
        X = fockobj.bket(sprintf('%dX',x));                                     %   x-coordinate ket
        Y = fockobj.bket(sprintf('%dY',y));                                     %   y-coordinate ket
        if     c==0; C = fockobj.bket('0C');                                    %   Cheese ket
        elseif c==1; C = fockobj.bket('1C');                                    %   ...
        else;        C = fockobj.bket('0C') + fockobj.bket('1C');               %   ...
        end                                                                     %   ...
        if     g==0; G = fockobj.bket('0G');                                    %   Water ket
        elseif g==1; G = fockobj.bket('1G');                                    %   ...
        else;        G = fockobj.bket('0G') + fockobj.bket('1G');               %   ...
        end                                                                     %   ...
        P  = [X Y C G] * [X Y C G]';                                            %   Make location projector
        obj.Oc(key) = P;                                                        %   Store in cache
      end                                                                       % <<
    end
    
    % TODO: Implement and write documentation comments
    function O=getVeridProj(obj)
      % Returns the veridicality operator.
      
      % TODO: Implement
    end

    % TODO: Write documentation comments
    function F=getFocusOp(obj,fea)
      % TODO: ...
      
      % Initialize                                                              % -------------------------------------
      fx = ~isemtpy(strfind(upper(fea),'X'));                                   % Focus on x-coordinate?
      fy = ~isemtpy(strfind(upper(fea),'Y'));                                   % Focus on y-coordinate?
      fc = ~isemtpy(strfind(upper(fea),'C'));                                   % Focus on cheese?
      fg = ~isemtpy(strfind(upper(fea),'G'));                                   % Focus on water?
      [xx,yy] = obj.mmp.getDim(); xx = xx+1; yy = yy+1;                         % Get coordianate ranges
      
      allX = 0; for i=1:xx; allX = allX + focjobj.bket(sprintf('%dX',i)); end
      allY = 0; for i=1:yy; allY = allY + focjobj.bket(sprintf('%dY',i)); end
      allC = 0; for i=0:1;  allC = allC + focjobj.bket(sprintf('%dC',i)); end
      allG = 0; for i=0:1;  allG = allG + focjobj.bket(sprintf('%dG',i)); end
      
      F = 0;
      for x=1:xx
        if fx; 
          xket = fockobj.bket(sprintf('%dX',x)); 
          xbra = xket';
        else
          xket = allX; 
          xbra = 1;
        end
        for y=1:yy
          if fy; 
            yket = fockobj.bket(sprintf('%dY',y)); 
            ybra = yket';
          else
            yket = allY; 
            ybra = 1;
          end
          for c=0:1
            if fc; 
              cket = fockobj.bket(sprintf('%dC',c)); 
              cbra = cket';
            else
              cket = allC; 
              cbra = 1;
            end
            for g=0:1
              if fg; 
                gket = fockobj.bket(sprintf('%dG',g)); 
                gbra = gket';
              else
                gket = allG; 
                gbra = 1;
              end
              F = F + [xket yket cket gket]*[xbra ybra cbra gbra];
              if ~fg; break; end
            end
            if ~fc; break; end
          end
          if ~fy; break; end
        end
        if ~fx; break; end
      end
    end

    % TODO: Write documentation comments
    function F=getUnfocusOp(obj,fea)
      F = obj.getFocusOp(fea);
      F = F';
    end

    function renderInnerStage(obj,lastAid)
    end
    
  end

  % == Protected constructors ==
  methods(Access=protected)

    function obj=Agent(type)
      % Creates a new agent
      %
      %   obj = Agent(type)
      %
      % arguments:
      %   type - The agent type, 'mouse' or 'cat'.
      %
      % returns:
      %   obj  - The mouse object.

      % Initialize mouse serial number                                          % -------------------------------------
      persistent serialNo;                                                      % Agent serial number
      if isempty(serialNo); serialNo = 0; end                                   % Initialize agent serial number
      serialNo = serialNo+1;                                                    % Get serial number of this agent
      obj.id = sprintf('%s#%d',type,serialNo);                                  % Create identifier of this agent

      % Initialize states, actions and operator cache                           % -------------------------------------
      obj.w   = 0;                                                              % Initial world state
      obj.phi = 0;                                                              % Initial mouse state
      obj.act = containers.Map;                                                 % Create actions map
      obj.Oc  = containers.Map;                                                 % Create operator cache
      
      % Initialize MouseMazePlot of inner stage                                 % -------------------------------------
      figure('Name',sprintf('Inner stage of %s',obj.id));                       % Create inner stage plot figure
      obj.mmp = MouseMazePlot([0 0]);                                           % Create inner stage plot

    end

  end
  
end

% EOF