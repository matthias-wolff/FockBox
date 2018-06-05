 classdef fock
  % Constants and static utility functions for the Fock space toolbox.
  %
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %   Peter beim Graben, BTU Cottbus-Senftenberg
  %
  % See also fockbasis, fockobj

  % == Settings ==
  properties(Constant)
    
    % Enables or disables additional checks (slower but safer).
    CHECK = true;
    
  end

  % == Constants ==
  properties(Constant)
      
    % Tensor product delimiter of base identifiers ('&#8855;').
    C_TDEL = char(8855);
    
    % Alternative tensor product delimiter of base identifiers ('&#167;').
    C_TDELA = char(167);
    
    % Identifier of vacuum base ('&#949;').
    C_BVAC = char(949);
    
    % Alternative identifier of vacuum base (ASCII,'').
    C_BVACA = '';
    
    % Row/column identidier delimiter ('&#10178;').
    C_RCDEL = char(10178);

  end

  % == Error message ids and messages ==
  properties(Constant)

    ERR = 'FOCK:err';
    
    % Error message indicating that an internal check has failed.
    ERR_CHECK = 'Internal check failed. %s';

    % Error message indicating that a function is not yet implemented.
    ERR_NOTIMPL = '%s not yet implemented.';
    
    % Error message indicating that an internal check has failed.
    ERR_MVERSION = 'Matlab version 9.4 (R2018a) or higher required for FockBox.';

    % Error message indicating a bad argument to a function.
    ERR_BADARG = 'Bad argument %s. %s';

    % Error message indicating a bad number of arguments to a function.
    ERR_BADNUMARGS = 'Bad number of arguments. %s';

    % Error message indicating too few arguments to a function.
    ERR_TOOMFEWARGS = 'Not enough arguments.';

    % Error message indicating too many arguments to a function.
    ERR_TOOMANYARGS = 'Too many arguments.';

    % Error message indicating a bad argument combination to a matrix product.
    ERR_BADMULT = 'Matrix multiplication %s not defined. %s';

    % Error message indicating a bad argument combination to a matrix product.
    ERR_BADINV = 'Inverse only defined for linear operators.';
    
    % Error message indicating that two sets are not disjoint.
    ERR_NOTDISJOINT = '%s must be disjoint';

    % Error message indicating that a base vector identifier was not found.
    ERR_BVEC_NOTFOUND = 'Base vector %s not found in %s';

    % Error message indicating that a <a href="matlab:helpPopup fockobj">fockobj</a> is invalid.
    ERR_OBJINVAL = "fockobj %s invalid. See "+fock.HREF('fockobj.getType')+".";

  end

  % == Fock space object types ==
  properties(Constant)

    % <a href="matlab:helpPopup fockobj">fockobj</a> is invalid.
    % See also fockobj.getType
    OBJ_INVALID = 0;

    % <a href="matlab:helpPopup fockobj">fockobj</a> is null (empty).
    % See also fockobj.getType
    OBJ_NULL = 1;

    % <a href="matlab:helpPopup fockobj">fockobj</a> is a scalar.
    % See also fockobj.getType
    OBJ_SCALAR = 2;

    % <a href="matlab:helpPopup fockobj">fockobj</a> is a ket.
    % See also fockobj.getType
    OBJ_KET = 3;

    % <a href="matlab:helpPopup fockobj">fockobj</a> is a bra.
    % See also fockobj.getType
    OBJ_BRA = 4;

    % <a href="matlab:helpPopup fockobj">fockobj</a> is a linear operator.
    % See also fockobj.getType
    OBJ_LOP = 5;

  end

  % == Public static methods
  methods(Static)

    function addFockBoxPaths()
      % Adds sub-folders in FockBox to Matlab path.
      %
      %   fock.addFockBoxPaths()
      %
      p = strrep(mfilename('fullpath'),mfilename,'');
      addpath([p 'example_utils']);
      addpath([p 'examples']);
      addpath([p 'external']);
    end
    
    function s=HREF(identifier)
      % Creates a hypertext link to a help page.
      % The returned string can be used for console printouts.
      % 
      %   s=fock.HREF(identifier)
      %
      % arguments:
      %     identifier - A function, class, or class member identifier.
      %
      % returns:
      %     s          - HTML text containing the hyperlink.
      s = sprintf(['<a href="matlab:helpPopup %s" style="font-weight:bold">'... % Create hypertext link
            '%s</a>'],identifier,identifier);                                   % ...
    end

    function C=vec2cell(V,type)
      % Converts a vector to a cell array of the specified type.
      %
      %   C=vec2cell(V,type)
      %
      % arguments:
      %   V    - The vector.
      %   type - The type, e.g. 'char', 'fockobj', etc.
      %
      % returns:
      %   C    - A cell array containing the elements of the vector or an
      %          empty cell array if the conversion failed.
      C={};                                                                     % Default return value 
      if ~isvector(V); return; end                                              % V is not a vector -> that's it
      if iscell(V)                                                              % V is a cell vector
        C = V;                                                                  %   Copy to result
      elseif ischar(V)                                                          % << V is a character vector >>
        C{1} = V;                                                               %   Store in a cell
      else                                                                      % << V is another vector >>
        C = cell(1,length(V));                                                  %   Pre-allocate resulting cell array
        for i=1:length(V); C{i}=V(i); end                                       %   Store elements of V as cells in C
      end                                                                       % <<
      if nargin>1                                                               % Do type checking/conversion
        assert(ischar(type)||(isstring(type)||isscalar(type)), ...              %   type must be a symbolic scalar
          fock.ERR_BADARG,'type','Must be a string or a character array');      %   ...
        type = string(type);                                                    %   Convert type to string
        for i=1:length(C)                                                       %   Loop over elements of C >>
          try                                                                   %     Try to convert cell to >>
            if type=="char"; C{i}=char(C{i}); end                               %       a char if type is "char"
            if type=="string"; C{i}=string(C{i}); end                           %       a string if type is "string"
          catch                                                                 %     << Type conversion failed >>
            C={}; return;                                                       %       Return empty vector
          end                                                                   %     <<
          if ~isa(C{i},type)                                                    %     Verify type of cell, if failed >>
            C={}; return;                                                       %       Return empty vector
          end                                                                   %     <<
        end                                                                     %   <<
      end                                                                       % <<
    end

    function [r,W]=plotVector(V,varargin)
      % Plots ket or bra vectors.
      % 
      %   [r,W]=plotVector(V)
      %   [r,W]=plotVector(___,Name,Value)
      %
      % arguments:
      %   V              - A vector of kets or bras.
      %   Name,Value     - Additional properties for plotting.
      %   ___            - Any other argument list.
      %
      % returns:
      %   r              - A struct of graphics objects that were drawn.
      %   W              - A matrix containing the (row) vectors actually
      %                    drawn, possibly after dimension reduction by PCA.
      %
      % additional properties:
      %   'mode'         - Drawing mode
      %     'arrows'     - Vector arrows from origin (default)
      %     'trajectory' - Vector trajectory
      %
      %   'axesLocation' - Location of coordinate axes
      %     'standard'   - Matlab standard (default)
      %     'origin'     - Axes through origin (experimental!)
      %
      %   'axesNames'    - Custom axes names (LaTeX format)
      %     {sx sy sz}   - Cell vector of names for x-, y-, and z-axis
      %
      %   'projection'   - Draw projection on coordinate plane (3D plots only)
      %     'off'        - Turn off (default)
      %     'xy'         - Draw projection on x-y plane
      %
      %   'projColor'    - Color of projection objects
      %     [r g b]      - An RGB array (no alpha!), default is [0 0.61 0.82]
      %
      %   'vecLabels'    - Vector labels (LaTeX format)
      %     {s1 ... sn}  - Cell vector with one label string per vector
      %
      %   'vecLabelOffs' - Position offsets of vector labels for fine-tuning
      %     {[dx1 dy1 dz1] ... [dxn dyn dzn]}
      %                  - Cell vector with one offset triplet per vector
      
      % Initialize                                                              % -------------------------------------
      V = fock.vec2cell(V,'fockobj');                                           % Convert vector set to cell array
      assert(~isempty(V),fock.ERR_BADARG,'V', ...                               % Check vector set V
        'Must be a scalar or a vector of type ''fockobj''.');                   % ...
      addpath(strrep(mfilename('fullpath'),mfilename,'external'));              % Add external folder to path
      
      % Customizable settings                                                   % -------------------------------------
      mode     = 0;             % Option 'mode'                                 % Vector array drawing mode (default)
      pmode    = 0;             % Option 'projection'                           % No projection on coord. plane (def.)
      caxNames = {};            % Option 'axesNames'                            % Custom axes names
      axOrigin = false;         % Option 'axesLocation'                         % Standard axes location (default)
      pcolor   = [0 0.61 0.82]; % Option 'projColor'                            % Projection color
      vLabels  = {};            % Option 'vecLabels'                            % Vector labels
      vlOffs   = {};            % Option 'vecLabelOffs'                         % Vector label position offsets
      for i=1:2:length(varargin)                                                % Loop over key-value pairs >>
        if strcmp(varargin{i},'mode')                                           %   Drawing mode >>
          if strcmp(varargin{i+1},'trajectory'); mode=1; end                    %     Trajectory
        end                                                                     %   <<
        if strcmp(varargin{i},'axesLocation')                                   %   Specification of axes locations >>
          axOrigin = strcmp(varargin{i+1},'origin');                            %     Through origin
        end                                                                     %   <<
        if strcmp(varargin{i},'axesNames')                                      %   Specification of axes names >>
          caxNames = varargin{i+1};                                             %     Get custom names
        end                                                                     %   <<
        if strcmp(varargin{i},'projection')                                     %   Projection on coordinate plane >>
          if strcmp(varargin{i+1},'xy'); pmode=1; end                           %     On x-y plane
        end                                                                     %   <<
        if strcmp(varargin{i},'projColor')                                      %   Color of projection objects >>
          pcolor = varargin{i+1};                                               %     Get color
        end                                                                     %   <<
        if strcmp(varargin{i},'vecLabels')                                      %   Specification of vector labels >>
          vLabels = varargin{i+1};                                              %     Get vector labels
        end                                                                     %   <<
        if strcmp(varargin{i},'vecLabelOffs')                                   %   Spec. of vec. label pos. offsets >>
          vlOffs = varargin{i+1};                                               %     Get position offsets
        end                                                                     %   <<
      end                                                                       % <<

      % Fix settings                                                            % -------------------------------------
      pcolor2 = pcolor+0.6*(1-pcolor);                                          % Light version of projection color
      tlN = {'Color', 'LineWidth', 'marker', 'MarkerSize'};                     % Trajectory line attributes
      tlV = {'black', 1.2, '.', 16};                                            % ...
      plN = {'Color', 'LineStyle'};                                             % Projection line attributes
      plV = {[pcolor 0.2], '--'};                                               % ...
      pmN = {'Color', 'marker', 'MarkerSize'};                                  % Projection line foot point attributes
      pmV = {pcolor2, '.', 12};                                                 % ...
      cpN = {'FaceColor', 'EdgeColor'};                                         % Coordinate plane attributes
      cpV = {[pcolor2 0.5], [pcolor2 0.5]};                                     % ...

      % Prepare data points                                                     % -------------------------------------
      H = fockbasis(V);                                                         % Get enclosing Fock subspace
      dim = H.getDim()-1;                                                       % Get space dimension w/o vacuum base
      pts = zeros(length(V),dim);                                               % Pre-allocate data points array
      for i=1:length(V)                                                         % Loop over Fock space objects >>
        fobj = V{i};                                                            %   Get current object
        fobjt = fobj.getType();                                                 %   Get type
        assert(fobjt==fock.OBJ_KET||fobjt==fock.OBJ_BRA,fock.ERR_BADARG, ...    %   Must be a ket or a bra
          sprintf("V{%d}",i),'Must be a ket or a bra');                         %   ...
        if fobjt==fock.OBJ_KET; fobj = fobj'; end                               %   Transpose ket -> bra
        bra = H.realize(fobj);                                                  %   Realize bra
        assert(size(bra,1)==1&&size(bra,2)==dim+1);                             %   Check dimensions of realization
        bra = bra(1,2:dim+1);                                                   %   Remove vacuum comp. (1st column)
        pts(i,:) = bra;                                                         %   Store in data points array
      end                                                                       % <<
      axNames{1}=['$|$' H.getBvecId(2) '$\rangle$'];                            % Set name of x axis
      axNames{2}=''; axNames{3}='';                                             % Set default names of y and z axes
      if (dim>1); axNames{2}=['$|$' H.getBvecId(3) '$\rangle$']; end            % Min. 2D -> set name of y axis
      if (dim>2); axNames{3}=['$|$' H.getBvecId(4) '$\rangle$']; end            % 3D -> set name of z axis
      axNames = strrep(axNames,fock.C_TDEL,'$\otimes$');                        % Make LaTeX \otimes
      if dim<=2                                                                 % Space one-dimensional >>
        pts = [pts zeros(length(V),3-dim)];                                     %   Add zeros as y-/z-coordinates
        dim = 2;                                                                %   Set 2D plot
      elseif dim>3                                                              % Space more than three-dimensional >>
        dim = 3;                                                                %   Set 3D plot
        [~, pts, ~] = pca(pts,'Centered',false);                                %   Do PCA
        if size(pts,2)<=2                                                       %   One or two principal components >>
          pts = [pts zeros(length(V),3-size(pts,2))];                           %     Add zeros as y-/z-coordinates
          dim = 2;                                                              %     Its going to be a 2D plot, too
        elseif size(pts,2)>3                                                    %   << More than three princ. comps. >>
          pts = pts(:,1:3);                                                     %     Keep the first three only
        end                                                                     %   <<
        axNames{1} = '$|PC_1\rangle$';                                          %   Rename x axis
        axNames{2} = '$|PC_2\rangle$';                                          %   Rename y axis
        axNames{3} = '$|PC_3\rangle$';                                          %   Rename z axis
      end                                                                       % <<
      if ~isempty(caxNames)                                                     % Have custom axes labels >>
        axNames(1:length(caxNames)) = caxNames;                                 %   Override automatic choice
      end                                                                       % <<
        
      % Prepare axes                                                            % -------------------------------------
      mn = min(0,min(pts,[],1));                                                % Minimum /
      mx = max(0,max(pts,[],1));                                                % Maximum
      dd = mn==mx;                                                              % Tune minima and maxima for plot
      mn = mn-0.5*dd;                                                           % - Make min. and max. unequal
      mx = mx+0.5*dd;                                                           %   ...
      dd = mx-mn;                                                               % - Enlarge axes dimensions by 18%
      mn = mn-0.09*dd;                                                          %   ...
      mx = mx+0.09*dd;                                                          %   ...

      % Plot axes                                                               % -------------------------------------
      view(dim);                                                                % Set plot dimension
      axis equal; if axOrigin; axis tight; end                                  % Axes settings
      set(gcf,'color','w');                                                     % Set plot window background to white
      ax = gca;                                                                 % Get current axes
      ax.Box = 'off';                                                           % Turn off the surrounding box
      if dim==2                                                                 % 2D plot >>
        axis([ mn(1) mx(1) mn(2) mx(2) ]);                                      %   Create axes with limits
        if axOrigin                                                             %   Axes through origin >>
          ax.XAxisLocation = 'origin';                                          %     Move x axis to origin
          ax.YAxisLocation = 'origin';                                          %     Move y axis to origin
          pos = get(gcf, 'Position');                                           %     Get physical display position
          m  = max(pos(3),pos(4));                                              %     Get max. physical display dim.
          dx = (mx(1)-mn(1))/pos(3)*m/50;                                       %     Width of arrow in x units
          dy = dx;%(mx(2)-mn(2))/pos(4)*m/50;                                   %     Height of arrow in y units
          line([mx(1)-dx mx(1) mx(1)-dx], [-dy/2 0 dy/2],'Color','black');      %     Draw x-axis arrow
          line([-dx/2 0 dx/2], [mx(2)-dy mx(2) mx(2)-dy],'Color','black');      %     Draw y-axis arrow
        end                                                                     %   <<
      else                                                                      % << 3D plot >>
        if axOrigin                                                             %   Axes through origins >>
          axlab = axis3atOrigin([ mn(1) mx(1) mn(2) mx(2) mn(3) mx(3) ]);       %     Create axes with limits
          r.xlabel = axlab.xlabel;                                              %     Store x axis label in result
          r.ylabel = axlab.ylabel;                                              %     Store x axis label in result
          r.zlabel = axlab.zlabel;                                              %     Store x axis label in result
          r.xlabel.String=axNames{1}; r.xlabel.Interpreter='latex';             %     Label x axis
          r.ylabel.String=axNames{2}; r.ylabel.Interpreter='latex';             %     Label y axis
          r.zlabel.String=axNames{3}; r.zlabel.Interpreter='latex';             %     Label z axis
        else                                                                    %   << Standard axes locations >>
          axis([ mn(1) mx(1) mn(2) mx(2) mn(3) mx(3) ]);                        %     Create axes with limits
        end                                                                     %   <<
      end                                                                       % <<
      if dim==2||~axOrigin                                                      % Using built-in axes >>
        r.xlabel = xlabel(axNames{1}); r.xlabel.Interpreter='latex';            %   Label x axis, and store in result
        r.ylabel = ylabel(axNames{2}); r.ylabel.Interpreter='latex';            %   Label y axis, and store in result
        r.zlabel = zlabel(axNames{3}); r.zlabel.Interpreter='latex';            %   Label z axis, and store in result
        if axOrigin                                                             %   Axes through origin >>
          r.xlabel.Position = [mx(1) -(mx(2)-mn(2))/50];                        %     Reposition label of x-axis
          r.xlabel.HorizontalAlignment = 'left';                                %     Adjust text alignment
          r.xlabel.VerticalAlignment = 'top';                                   %     ...
          r.ylabel.Position = [-(mx(1)-mn(1))/50 mx(2)];                        %     Reposition label of y-axis
          r.ylabel.HorizontalAlignment = 'right';                               %     Adjust text alignment 
          r.ylabel.VerticalAlignment = 'middle';                                %     ...
          % Note: No need for z-label as 3D plots are handled by axis3atOrigin  %
        end                                                                     %   <<
      end                                                                       % <<
      
      % Plot data points                                                        % -------------------------------------
      if mode==1&&length(V)==1                                                  % Trajectory mode but only one point >>
        warning('Only one vector given. Falling back to arrow drawing mode.');  %   Warning
        mode = 0;                                                               %   Set vector arrow mode
      end                                                                       % <<
      if pmode>0&&dim==2                                                        % Projection selected but plot 2D >>
        warning('Plot is 2D. Turing off projection on coordinate plane.');      %   Warning
        pmode = 0;                                                              %   Set vector arrow mode
      end                                                                       % <<
      hold on;                                                                  % Hold plot on current axes
      if pmode==1                                                               % Project on x-y plane >>
        r.xyplane = rectangle('Position', ...                                   %   Draw x-y plane
          [mn(1)+0.09*dd(1) mn(2)+0.09*dd(2) mx(1)-mn(1)-0.18*dd(1) ...         %   ...        
           mx(2)-mn(2)-0.18*dd(2)], cpN, cpV);                                  %   ...
        for i=1:length(V)                                                       %   Loop over vectors >>
          if abs(pts(i,1))<dd(1)/100&&abs(pts(i,2))<dd(2)/100; continue; end    %     Do not print over z axis
          r.vec{i}.pline = line([pts(i,1) pts(i,1)], [pts(i,2) pts(i,2)], ...   %     Draw projection line on x-y plane
            [0 pts(i,3)], plN, plV);                                            %     ...
          r.vec{i}.pfoot = line([pts(i,1) pts(i,1)], [pts(i,2) pts(i,2)], ...   %     Draw projection line foot point
            [0 0], pmN, pmV);                                                   %     ...
        end                                                                     %   <<
      end                                                                       % <<
      if mode==0                                                                % Vector arrow drawing mode >>
        for i=1:length(V)                                                       %   Loop over vectors >>
          r.vec{i}.arrow = arrow([0 0 0],[pts(i,1) pts(i,2) pts(i,3)], ...      %     Plot vector arrow from origin
            'BaseAngle',60,'LineWidth',1.2);                                    %     ...
        end                                                                     %   <<
      else                                                                      % << Trajectory drawing mode >>
        for i=2:length(V)                                                       %   Loop over vectors >>
          r.vec{i}.line = line([pts(i-1,1) pts(i,1)], [pts(i-1,2) pts(i,2)], ...%     Draw line
            [pts(i-1,3) pts(i,3)], tlN, tlV);                                   %     ...
        end                                                                     %   <<
      end                                                                       % <<
      
      % Plot data point labels                                                  % -------------------------------------
      for i=1:length(V)                                                         % Loop over vectors >>
        if i<=length(vLabels)                                                   %   Custom vector labels >>
          s = char(vLabels{i});                                                 %     Use them
        else                                                                    %   << No custom vector labels >>
          s = char(V{i}.name);                                                  %     Use name of fockobj
        end                                                                     %   <<
        if strcmp(s,""); s=''; end                                              %   Unnamed -> make empty name
        l = sqrt(sum(pts(i,:).^2));                                             %   Compute length of vector
        pos = pts(i,:).*(1+0.09*dd./l);                                         %   Compute label position
        offs = [0 0 0];                                                         %   Label offsets dx, dy, and zz
        if i<=length(vlOffs)                                                    %   Have custom offset(s) >>
          offs(1,1:length(vlOffs{i})) = vlOffs{i};                              %     Apply
        end                                                                     %   <<
        pos = pos+offs;                                                         %   Offset vector labels
        r.vec{i}.label = text(pos(1), pos(2), pos(3), s, ...                    %   Draw label (even if empty!)
          'HorizontalAlignment','center','VerticalAlignment','middle');         %   ...
        if any(strfind(s,'$')); r.vec{i}.label.Interpreter = 'latex'; end       %   HACK: auto-detect LaTeX mode
      end                                                                       % <<

      % Finish plot                                                             % -------------------------------------
      hold off;                                                                 % Release hold
      W = pts;                                                                  % Return actually drawn data points

    end

    function setVectorColor(v,c)
      % Sets the color of a vector arrow plotted by plotVector.
      %
      %   setArrowColor(v,c)
      %
      % arguments:
      %   v - A vector arrow plotted by <a href="matlab:helpPopup plotVector">plotVector</a>.
      %   c - A color as defined in <a href="matlab:helpPopup patch">patch</a>.
      v.arrow.FaceColor=c; v.arrow.EdgeColor=c;
      v.label.Color=c;
    end
    
  end

end

% EOF