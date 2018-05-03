classdef fockpst < fockobj
  % A phrase structure tree in Chomsky normal form (subclass of fockobj).
  %
  % example:
  %   Phrase structure tree for "Colorless green ideas sleep furiously":
  %
  %   obj = fockpst(              ...
  %     { 'S'                     ...
  %       { 'NP'                  ...
  %         { 'A' 'colorless' }   ...
  %         { 'NP'                ...
  %           { 'A' 'green' }     ...
  %           { 'N' 'ideas' }     ...
  %         }                     ...
  %       }                       ...
  %       { 'VP'                  ...
  %         { 'V' 'sleep' }       ...
  %         { 'Adv' 'furiously' } ...
  %       }                       ...
  %     });
  %
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %   Peter beim Graben, BTU Cottbus-Senftenberg
  %
  % See also fockobj  

  properties(Constant)
    
    % The mother role base.
    R_M  = fockobj.bket('^');
    
    % The left-child role base.
    R_LC = fockobj.bket('/');
    
    % The right-child role base.
    R_RC = fockobj.bket('\');
    
  end
  
  methods
    
    function obj = fockpst(ca)
      % Creates a ket representing a phrase structure tree.
      %
      %   obj=fockpst(ca)
      %
      % arguments:
      %   ca  - A cell array defining a phrase structre tree in Chomsky
      %         normal form.
      %
      % returns:
      %   obj - A ket representing the phrase structure tree.
      f = fockpst.fromCellArray(ca);                                           % Create fockobj
      obj = obj@fockobj();                                                      % Invoke superclass constructor
      obj.data = f.data;                                                        % Copy data
    end
    
  end
  
  methods(Static, Access=private)
   
    % TODO: Write doc and line comments!
    function obj=fromCellArray(ca)
      % TODO: ...
      
      if isempty(ca)
        obj = fockpst.R_M;
        return;
      end
      
      assert(iscell(ca)&&isvector(ca), ...
        fock.ERR_BADARG,'ca','Must be a cell vector');
      assert(length(ca)==2||length(ca)==3, ...
        fock.ERR_BADARG,'ca','Sub-cells must be vectors with 2 or 3 elements');
      assert(ischar(ca{1})||isstring(ca{1}), ...
        fock.ERR_BADARG,'ca','First element of sub-cells must be a string');
      
      obj = [ fockobj.bket(ca{1}) fockpst.R_M ];
      if iscell(ca{2})
        obj = obj + [ fockpst.fromCellArray(ca{2}) fockpst.R_LC ];
      else
        obj = obj + [ fockobj.bket(ca{2}) fockpst.R_LC ];
      end
      if (length(ca)==3)
        if iscell(ca{3})
          obj = obj + [ fockpst.fromCellArray(ca{3}) fockpst.R_RC ];
        else
          obj = obj + [ fockobj.bket(ca{3}) fockpst.R_RC ];
        end
      end
    end
    
  end
  
end

% EOF