classdef Agent < handle
  % An agent living in a maze world.
  %
  % See also: Maze
  
  %% == Properties ==
  properties(SetAccess=protected)
    
    % Unique identifier of agent.
    id;
    
    % MouseMazePlot visualizing the agent's inner stage.
    mmp;

    % Flag indicating that inner stage plot needs to be redrawn.
    mmpStale;
    
    % Flag indicating that redrawing inner stage plot is deferred.
    mmpDefer;
    
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
  
  %% == Instance getters ==
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
  
  %% == Setters and getters ==
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
    
    function deferPlot(obj,defer)
      % Toggles deferring update of the inner stage plot.
      %
      %    obj.deferPlot(defer)
      %
      % arguments:
      %    obj   - The agent.
      %    defer - True to defer updating the plot, false to enable update.
      %
      % remarks:
      %    Invoking deferPlot(false) does not trigger an update of the plot. 
      %    The plot will be updated by the next invocation of the 
      %    innerStageChanged(...) function.
      %
      % See also: innerStageChanged
      
      obj.mmpDefer = defer~=0;
    end
    
  end
  
  %% == Perception-action cycle
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
      
      a.name   = sprintf('a_%s',aid);                                           % Make world action name
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
      changed = 0;                                                              % Inner stage changed by perception
      
      % Update trial action operator                                            % -------------------------------------
      if nargin>=3                                                              % Action caused observation >>
        action = obj.act(aid);                                                  %   Fetch action
        Oold = action.O;                                                        %   Remember old trial action
        action.O = action.O + (1-o'*action.O*obj.phi)*(o*obj.phi');             %   Update trial action
        action.O.name = sprintf('O_%s',aid);                                    %   Set trial action name
        obj.act(aid) = action;                                                  %   Store action
        if Oold~=action.O; changed = 1; end                                     %   Inner stage changed
        disp(action.O);                                                         %   Console log
      else                                                                      % << Other observation >>
        aid = '';                                                               %   aid <- no-operation
      end                                                                       % <<

      % Update agent and world states                                           % -------------------------------------
      phiold  = obj.phi;                                                        % Remember old agent state
      wold    = obj.w;                                                          % Remember old world state
      obj.phi = o;                                                              % Update agent state
      obj.w   = obj.w + (1-obj.w'*o)*o;                                         % Update world state
      if phiold~=obj.phi || wold~=obj.w; changed = 1; end                       % Inner stage changed

      % Aftermath                                                               % -------------------------------------
      if changed                                                                % Inner stage changed by perception >>
        if obj.mmpStale==0                                                      %   This was the only change >>
          [x,y,p] = obj.locate();                                               %     Locate agent
          if p==1                                                               %       Agent loc. fully determined >>
            obj.innerStageChanged(aid,x,y);                                     %       Incremental update of plot
          else                                                                  %     << Otherwise >>
            obj.innerStageChanged(aid);                                         %       Full update of plot
          end                                                                   %     <<
        else                                                                    %   << There were more changes >>
          obj.innerStageChanged(aid);                                           %     Full update of plot
        end                                                                     %   <<
      end                                                                       % <<
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

    function sexplore(obj,rx,ry,fast)
      % Performs a systematic exploration.
      %
      %   obj.sexplore(rx,ry)
      %
      % arguments:
      %   obj  - The agent.
      %   rx   - Range of x-coordinates to explore, e.g. 1:10.
      %   ry   - Range of y-coordinates to explore.
      %   fast - true for speed-up (skips most of the drawing, optional)
      
      if nargin<4; fast = 0; end                                                % Default fast option (off)
      aids = obj.act.keys;                                                      % Get action ids
      for y=ry                                                                  % Loop over y-coordinates >>
        for x=rx                                                                %   Loop over x-coordinates >>
          if fast; obj.deferPlot(1); end                                        %     Fast mode: defer drawing
          for i=1:length(aids)                                                  %     Loop over actions >>
            obj.teleport(x,y);                                                  %       Teleport agent to coordinates
            obj.action(aids{i});                                                %       Perform action
          end                                                                   %     <<
          if fast                                                               %     Fast mode >>
            obj.deferPlot(0);                                                   %       Set defer drawing off
            obj.innerStageChanged('',x,y);                                      %       Draw inner stage
            if obj.isMouse(); obj.mmp.moveMouse(x,y,200); end                   %       Draw mouse agent
            % TODO: Other agents
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      
    end

  end
  
  %% == Inner stage ==
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
      
      P = obj.getLocationProj(x,y);                                             % Get location projector
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

    function P=getLocationProj(obj,x,y,c,g)
      % Returns a location projector.
      %
      %   P = obj.getLocationProj(x,y,c,g)
      %
      % arguments:
      %   obj - The agent.
      %   x   - The x-coordinate of the location.
      %   y   - The y-coordinate of the location.
      %   c   - Cheese flag: 0 - no cheese, 1 - cheese, otherwise - dont' care (optional)
      %   g   - Water flag: 0 - no water, 1 - water, otherwise - dont' care (optional)
      %
      % returns:
      %   P   - The location projector.

      % Initialize                                                              % -------------------------------------
      if nargin<4; c=-1; end; if c~=0&&c~=1; c=-1; end                          % Default cheese flag, normalize
      if nargin<5; g=-1; end; if g~=0&&g~=1; g=-1; end                          % Default water flag, normalize

      % Cache retrieval                                                         % -------------------------------------
      key = sprintf('PL(%d,%d,%d,%d)',x,y,c,g);                                 % Make operator cache key
      if obj.Oc.isKey(key)                                                      % Veridicality projector in cache >>
        P = obj.Oc(key);                                                        %   Retrieve veridicality projector
        return;                                                                 %   That's it...
      end                                                                       % <<

      % Build operator                                                          % -------------------------------------
      X  = fockobj.bket(sprintf('%dX',x));                                      % x-coordinate ket
      Y  = fockobj.bket(sprintf('%dY',y));                                      % y-coordinate ket
      C0 = fockobj.bket('0C');                                                  % No-cheese ket
      C1 = fockobj.bket('1C');                                                  % Cheese ket
      G0 = fockobj.bket('0G');                                                  % No-water ket
      G1 = fockobj.bket('1G');                                                  % Water ket
      P = 0;                                                                    % Initial projector
      if (c==0||c==-1)&&(g==0||g==-1); P = P + [X Y C0 G0] * [X Y C0 G0]'; end  % No cheese and no water
      if (c==1||c==-1)&&(g==0||g==-1); P = P + [X Y C1 G0] * [X Y C1 G0]'; end  % Cheese but no no water
      if (c==0||c==-1)&&(g==1||g==-1); P = P + [X Y C0 G1] * [X Y C0 G1]'; end  % No cheese but water
      if (c==1||c==-1)&&(g==1||g==-1); P = P + [X Y C1 G1] * [X Y C1 G1]'; end  % Cheese and water
      P.name = key;                                                             % Name projector

      % Cache location projector                                                % -------------------------------------
      obj.Oc(key) = P;                                                          % Write to operator cache
    end

    function P=getVeridicalityProj(obj)
      % Returns the veridicality projector.
      % 
      %   P = obj.getVeridicalityProj()
      %
      % arguments:
      %   obj - The agent.
      %
      % returns:
      %   P   - The veridicality projector.

      % Cache retrieval                                                         % -------------------------------------
      if obj.Oc.isKey('PV')                                                     % Veridicality projector in cache >>
        P = obj.Oc('PV');                                                       %   Retrieve veridicality projector
        return;                                                                 %   That's it...
      end                                                                       % <<

      % Build operator                                                          % -------------------------------------
      [xx,yy] = obj.mmp.getDim();                                               % Get estimated coordinate ranges
      wBra = obj.w';                                                            % Get <w|
      P = 0;                                                                    % Initialize projector
      for x=1:xx                                                                % Loop over x-coordinates >>
        X = fockobj.bket(sprintf('%dX',x));                                     %   Make |xX>
        for y=1:yy                                                              %   Loop over y-coordinates >>
          Y = fockobj.bket(sprintf('%dY',y));                                   %     Make |yY>
          for c=0:1                                                             %     Loop over cheese values >>
            C = fockobj.bket(sprintf('%dC',c));                                 %       Make |cC>
            for g=0:1                                                           %       Loop over water values
              G = fockobj.bket(sprintf('%dG',g));                               %         Make |gG>
              B = [X Y C G];                                                    %         Make "basic situation" ket
              P = P + (wBra*B)*(B*B');                                          %         Proj. on |w> and add
            end                                                                 %       <<
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      P.name = 'PV';                                                            % Name projector

      % Cache veridicality projector                                            % -------------------------------------
      obj.Oc('PV') = P;                                                         % Write to operator cache

    end

    function F=getFocusOp(obj,fea)
      % Returns a focus operator.
      %
      %   F = obj.getFocusOp(fea)
      %
      % arguments:
      %   obj - The agent.
      %   fea - Any combination of
      %         * 'X': focus on x-coordinate feature,
      %         * 'Y': focus on y-coordinate feature,
      %         * 'C': focus on cheese feature, and
      %         * 'G': focus on water feature,
      %         but not the empty char array ''.
      %
      % returns:
      %   F   - The focus operator.
      %
      % See also: getUnfocusOp

      % Initialize                                                              % -------------------------------------
      key = '';                                                                 % Normalized version of fea
      fx = contains(upper(fea),'X'); if fx; key = [key 'X']; end                % Focus on x-coordinate?
      fy = contains(upper(fea),'Y'); if fy; key = [key 'Y']; end                % Focus on y-coordinate?
      fc = contains(upper(fea),'C'); if fc; key = [key 'C']; end                % Focus on cheese?
      fg = contains(upper(fea),'G'); if fg; key = [key 'G']; end                % Focus on water?
      key = sprintf('F(%s)',key);                                               % Make operator cache key

      % Cache retrieval                                                         % -------------------------------------
      if obj.Oc.isKey(key)                                                      % Focus operator in cache >>
        F = obj.Oc(key);                                                        %   Retrieve focus operator
        return;                                                                 %   That's it...
      end                                                                       % <<

      % Build operator                                                          % -------------------------------------
      [xx,yy] = obj.mmp.getDim();                                               % Get estimated coordinate ranges
      allX = 0; for i=1:xx; allX = allX + fockobj.bket(sprintf('%dX',i)); end   % Ket representing all x-coordinates
      allY = 0; for i=1:yy; allY = allY + fockobj.bket(sprintf('%dY',i)); end   % Ket representing all y-coordinates
      allC = 0; for i=0:1;  allC = allC + fockobj.bket(sprintf('%dC',i)); end   % Ket representing cheese and no cheese
      allG = 0; for i=0:1;  allG = allG + fockobj.bket(sprintf('%dG',i)); end   % Ket representing water and no water
      F = 0;                                                                    % Initial focus operator
      for x=1:xx                                                                % Loop over x-coordinates >>
        if fx                                                                   %   Focus on x >>
          xket = fockobj.bket(sprintf('%dX',x));                                %     Ket factor <- |xX>
          xbra = xket';                                                         %     Bra factor <- <xX|
        else                                                                    %   << Do not focus on x >>
          xket = 1;                                                             %     Ket factor <- 1
          xbra = allX';                                                         %     Bra factor <- sum_i <iX|
        end                                                                     %   <<
        for y=1:yy                                                              %   Loop over y-coordinates >>
          if fy                                                                 %     Focus on y >>
            yket = fockobj.bket(sprintf('%dY',y));                              %       Ket factor <- |yY>
            ybra = yket';                                                       %       Bra factor <- <yY|
          else                                                                  %     << Do not focus on y >>
            yket = 1;                                                           %       Ket factor <- 1
            ybra = allY';                                                       %       Bra factor <- sum_i <iY|
          end                                                                   %     <<
          for c=0:1                                                             %     Loop over cheese values >>
            if fc                                                               %       Focus on cheese >>
              cket = fockobj.bket(sprintf('%dC',c));                            %         Ket factor <- |cC>
              cbra = cket';                                                     %         Bra factor <- <cC|
            else                                                                %       << Do not focus on cheese >>
              cket = 1;                                                         %         Ket factor <- 1
              cbra = allC';                                                     %         Bra factor <- sum_i <iC|
            end                                                                 %       <<
            for g=0:1                                                           %       Loop over water values >>
              if fg                                                             %         Focus on water >>
                gket = fockobj.bket(sprintf('%dG',g));                          %           Ket factor <- |gG>
                gbra = gket';                                                   %           Bra factor <- <gG|
              else                                                              %         << Do not focus on water >>
                gket = 1;                                                       %           Ket factor <- 1
                gbra = allG';                                                   %           Bra factor <- sum_i <iG|
              end                                                               %         <<
              F = F + [xket yket cket gket]*[xbra ybra cbra gbra];              %         Aggregate focus operator
              if ~fg; break; end                                                %         No more water summands
            end                                                                 %       <<
            if ~fc; break; end                                                  %       No more cheese summands
          end                                                                   %     <<
          if ~fy; break; end                                                    %     No more y summands
        end                                                                     %   <<
        if ~fx; break; end                                                      %   No more x summands
      end                                                                       % <<
      F.name = key;                                                             % Name operator

      % Cache focus operator                                                    % -------------------------------------
      obj.Oc(key) = F;                                                          % Write to operator cache

    end

    function F=getUnfocusOp(obj,fea)
      % Returns a generic unfocus operator. Multiply veridicality projector
      % with return value to obtain the respective veridical unfocus operator.
      %
      %   F = obj.getUnfocusOp(fea)
      %
      % arguments:
      %   obj - The agent.
      %   fea - Any combination of
      %         * 'X': focus on x-coordinates,
      %         * 'Y': focus on y-coordinates,
      %         * 'C': focus on cheese, and
      %         * 'G': focus on water,
      %         but not the empty char array ''.
      %
      % returns:
      %   F   - The unfocus operator.
      %
      % See also: getFocusOp, getVeridicalityProj

      F = obj.getFocusOp(fea);                                                  % Get respective focus operator
      name = F.name;                                                            % Get name of focus operator
      F = F';                                                                   % Return conjugate transpose
      F.name = strrep(name,'F(','Ug(');                                         % Name unfocus operator
    end

    function innerStageChanged(obj,aid,x,y)
      % Invoked when the inner stage has changed. Removes invalidated 
      % operators from the operator cache and redraws the plot.
      %
      %    obj.innerStageChanged(aid,x,y)
      %    obj.innerStageChanged(aid)
      %    obj.innerStageChanged()
      %
      % arguments:
      %    obj - The agent.
      %    aid - Action that caused the change, '' for no action (optional).
      %    x   - x-coordinate for incremental update of inner stage plot (optional).
      %    y   - y-coordinate for incremental update of inner stage plot (optional).
      %
      % See also: deferPlot
      
      % Initialize                                                              % -------------------------------------
      if nargin<2; aid=''; end                                                  % Default action id: '' (no-operation)
      
      % Remove invalidated operators from cache                                 % -------------------------------------
      if obj.Oc.isKey('PV'); obj.Oc.remove('PV'); end                           % Remove veridicality proj. from cache
      keys = obj.Oc.keys(); obj.Oc.remove(keys(startsWith(keys,'OF(')));        % Remove focus operators from cache
      
      % Update plot                                                             % -------------------------------------
      if obj.mmpDefer                                                           % Plotting deferred >>
        obj.mmpStale = obj.mmpStale + 1;                                        %   Increment stale counter
        return;                                                                 %   Slip plot update
      end                                                                       % <<
      
      % - Initialize                                                            % - - - - - - - - - - - - - - - - - - -
      [xx,yy] = obj.mmp.getDim();                                               % Get (known) maze dimension
      fprintf('\nEstimated maze size: %d x %d\n',xx,yy);                        % Console log
      if nargin>=4                                                              % Incremental update >>
        if obj.isMouse(); obj.mmp.removeAllMice(); end                          %   Remove mouse agent from plot
        % TODO: Other agents
        rx = x:x; ry=y:y;                                                       %   Coordinate range to update
      else                                                                      % << Full update >>
        obj.mmp.clear();                                                        %   Clear plot
        rx = 1:xx+1; ry=1:yy+1;                                                 %   Coordinate range to update
      end                                                                       % <<
      aids = ['N', 'E', 'S', 'W'];                                              % Movement action ids
      h = 200;                                                                  % Default agent heading
      if     aid=='N'; h = 0;                                                   % Agent heading for north action
      elseif aid=='E'; h = 90;                                                  % Agent heading for east action
      elseif aid=='S'; h = 180;                                                 % Agent heading for south action
      elseif aid=='W'; h = -90;                                                 % Agent heading for west action
      end                                                                       %

      % - Update                                                                % - - - - - - - - - - - - - - - - - - -
      for y=ry                                                                  % Loop over y-coordinates >>
        for x=rx                                                                %   Loop over x-coordinates >>
          PL = obj.getLocationProj(x,y);                                        %     Get location projector (x,y)
          pl = obj.w'*PL*obj.w;                                                 %     Probability of location (x,y)
          if pl==0; continue; end                                               %     Skip unknown location
          obj.mmp.addFloorTile(x,y,pl);                                         %     Add place
          pc = obj.w'*obj.getLocationProj(x,y,1)*obj.w;                         %     Prob. of cheese at location (x,y)
          obj.mmp.addCheese(x,y,pc);                                            %     Add cheese
          pg = obj.w'*obj.getLocationProj(x,y,-1,1)*obj.w;                      %     Prob. of water at location (x,y)
          obj.mmp.addWater(x,y,pg);                                             %     Add water
          pa = obj.phi'*PL*obj.phi;                                             %     Prob. of agent at location (x,y)
          if pa>0                                                               %     Agent probability >0 >>
            if obj.isMouse(); obj.mmp.addMouse(x,y,pa,h); end                   %       Add mouse agent
            % TODO: Other agents
          end                                                                   %     <<
          X = fockobj.bket(sprintf('%dX',x));                                   %     x-coordinate ket
          Y = fockobj.bket(sprintf('%dY',y));                                   %     y-coordinate ket
          C = (1-pc)*fockobj.bket('0C') + pc*fockobj.bket('1C');                %     Cheese ket at (x,y)
          G = (1-pg)*fockobj.bket('0G') + pg*fockobj.bket('1G');                %     Water ket at (x,y)
          phixy = [X Y C G];                                                    %     Agent ket at (x,y)
          for i=1:length(aids)                                                  %     Loop over movement actions >>
            aid = aids(i);                                                      %       Get action identifier
            if ~obj.act.isKey(aid); continue; end                               %       Skip unknown action ids
            action = obj.act(aid);                                              %       Get action
            if phixy.eq(action.O*phixy)                                         %       Trial action ineffective >>
              obj.mmp.addWall(x,y,aid,0.3);                                     %         Add wall
            end                                                                 %       <<
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
      
      % - Redraw                                                                % - - - - - - - - - - - - - - - - - - -
      drawnow;                                                                  % Update plots  
      obj.mmpStale = 0;                                                         % Reset stale counter
    end
    
  end

  %% == Protected constructors ==
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
      obj.mmpStale = 0;                                                         % Plot up to date
      obj.mmpDefer = 0;                                                         % Immediately redraw plot

    end

  end
  
end

% EOF