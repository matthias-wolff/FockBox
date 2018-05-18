classdef Mouse < handle
  %MOUSE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(SetAccess=protected)
    
    % Unique identifier of agent.
    id;
    
    % MouseMazePlot visualizing the mouse's inner stage.
    mmp;

    % The world state
    w;
    
    % The mouse state
    phi;
    
    % The action repertiore
    act;

  end
  
  properties(SetAccess={?Maze})
    
    % The Maze world agent is living in.
    world;
    
  end
  
  methods

    function obj=Mouse()
      % Creates a new mouse agent
      %
      %   obj = Mouse()
      %
      % returns:
      %   obj   - The mouse object.

      % Initialize mouse serial number                                          % -------------------------------------
      persistent serialNo;                                                      % Mice serial number
      if isempty(serialNo); serialNo = 0; end                                   % Initialize mice serial number
      serialNo = serialNo+1;                                                    % Get serial number of this mouse
      obj.id = sprintf('mouse%d',serialNo);                                     % Store serial number of this mouse

      % Initialize states and actions                                           % -------------------------------------
      obj.w   = 0;                                                              % Initial world state
      obj.phi = 0;                                                              % Initial mouse state
      obj.act = containers.Map;                                                 % Create actions map
      
      % Initialize MouseMazePlot of inner stage                                 % -------------------------------------
      figure('Name',obj.id);                                                    % Create inner stage plot figure
      obj.mmp = MouseMazePlot([0 0]);                                           % Create inner stage plot

    end

    % TODO: Write documentation comments
    function addAction(obj,aid,a)
      % TODO: ...
      
      action.a = a;                                                             % Create world action
      action.O = 0;                                                             % Create trial action
      obj.act(aid) = action;                                                    % Add action
    end

    % TODO: Write documentation comments
    function o=teleport(obj,x,y)
      % TODO: ...

      fprintf('__________________________________________________________\n');
      fprintf('TELEPORT agent %s to (%d,%d):\n',obj.id,x,y);
      disp(a);

      o = perceive(obj.world.teleport(obj,x,y));                                % Teleport and preceive result

      fprintf('\n-> Agent state:\n'); disp(obj.phi);
      fprintf('\n-> World state:\n'); disp(obj.w);
      fprintf('__________________________________________________________\n');
    end
    
    % TODO: Write documentation comments
    function o=action(obj,aid)
      % TODO: ...

      a = obj.act(aid).a;                                                       % Get action semantics
      fprintf('___________________________________________________________\n'); % Console log
      fprintf('ACTION "%s" of agent %s:\n',aid,obj.id);                         % Console log
      disp(a);                                                                  % Console log
      o = obj.perceive(obj.world.action(obj,a));                                % Execute action and preceive result
      fprintf('\n> Agent state:\n'); disp(obj.phi);                             % Console log
      fprintf('\n> World state:\n'); disp(obj.w);                               % Console log
      fprintf('___________________________________________________________\n'); % Console log
    end

    % TODO: Write documentation comments
    function o=perceive(obj,o)
      % TODO: ...

      obj.phi = o;                                                              % Update mouse state
      obj.w = obj.w + (1-obj.w'*o)*o;                                           % Update maze state
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
      % Perform systematic exploration.
      
      % TODO: Implement
    end
    
    % TODO: Implement and write documentation comments
    function O=getVeridicalOp(obj)
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

end

% EOF