classdef fockfvr < fockobj
  % A feature-value relation (subclass of fockobj).
  %
  % example:
  %   fvr = fockfvr(   ...
  %     { 'Maze'       ...
  %       { 'Location' ...
  %         { 'X'      ...
  %           { 1 }    ...
  %         }          ... 
  %         { 'Y'      ...
  %           { 1 }    ...
  %         }          ...
  %       }            ...
  %       { 'Cheese'   ...
  %         { 1 }      ...
  %       }            ...
  %     });  
  % 
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %   Peter beim Graben, BTU Cottbus-Senftenberg
  %
  % See also fockobj  

  methods
    
    function obj = fockfvr(ca)
      % Creates a ket representing a feature-value relation.
      %
      %   obj=fockpst(ca)
      %
      % arguments:
      %   ca  - A cell array defining a feature-value relation.
      %
      % returns:
      %   obj - A ket representing the feature-value relation.
      f = fockfvr.fromCellArray(ca);                                            % Create fockobj
      obj = obj@fockobj();                                                      % Invoke superclass constructor
      obj.data = f.data;                                                        % Copy data
    end
        
  end
  
  methods(Static, Access=private)
   
    % TODO: Write doc and line comments!
    function obj=fromCellArray(ca)
      % TODO: ...
      
      if isempty(ca)
        obj = 0;
        return;
      end
      
      assert(~iscell(ca)||isvector(ca), ...
        fock.ERR_BADARG,'ca','Must be a string or a cell vector');
      assert(~isempty(ca), ...
        fock.ERR_BADARG,'ca','Sub-cells must not be empty');

      nodeid = '';
      subfvr = [];
      for i=1:length(ca)
        if ~iscell(ca{i})
          nodeid = [ nodeid char(string(ca{i})) ];
        else
          subfvr = [ subfvr; fockfvr.fromCellArray(ca{i}) ];
        end
      end
      if isempty(nodeid)
        warning('Missing node identifier. Assuming ''_''.');
        nodeid = '_';
      end
      node = fockobj.bket(nodeid);
      if isempty(subfvr)
        obj= node;
      else
        obj = [ node subfvr(1) ];
        for i=2:length(subfvr)
          obj = obj + [ node subfvr(i) ];
        end
      end
      
    end
    
  end

end

% EOF