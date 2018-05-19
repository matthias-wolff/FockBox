classdef Mouse < handle
  %MOUSE Summary of this class goes here
  %   Detailed explanation goes here
  
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

    function p=isat(obj,x,y)
    end

    function [x y p]=locate(obj)
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

      if nargin>=3
        % TODO: Compare o and obj.phi -> update trial action operator
      end
      
      obj.phi = o;                                                              % Update mouse state
      obj.w = obj.w + (1-obj.w'*o)*o;                                           % Update maze state
      fprintf('\n> Agent state:\n'); disp(obj.phi);                             % Console log
      fprintf('\n> World state:\n'); disp(obj.w);                               % Console log
    end

    % TODO: Implement and write documentation comments
    function rexplore(obj,n)
      % Perform random exploration.
      
      aids = obj.act.keys;                                                      % Get action ids
      for i=1:n                                                                 % Random exploration iterations >>
        aid = aids{randi(length(aids))};                                        %   Randomly select action id
        obj.action(aid);                                                        %   Perform action
      end                                                                       % <<
      
    end

    % TODO: Implement and write documentation comments
    function sexplore(obj,rx,ry)
      % Performs a systematic exploration.
      
      aids = obj.act.keys;                                                      % Get action ids
      for x=rx                                                                  % Loop over x-coordinates >>
        for y=ry                                                                %   Loop over y-coordinates >>
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

      % Initialize states and actions                                           % -------------------------------------
      obj.w   = 0;                                                              % Initial world state
      obj.phi = 0;                                                              % Initial mouse state
      obj.act = containers.Map;                                                 % Create actions map
      
      % Initialize MouseMazePlot of inner stage                                 % -------------------------------------
      figure('Name',obj.id);                                                    % Create inner stage plot figure
      obj.mmp = MouseMazePlot([0 0]);                                           % Create inner stage plot

    end

  end
  
end

% EOF