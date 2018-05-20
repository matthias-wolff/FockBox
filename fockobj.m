classdef fockobj < matlab.mixin.CustomDisplay
  % A Fock space object is scalar, a ket, a bra, or a linear operator living in 
  % a Fock space. Kets are created by the static 'bket' method
  %
  %   X = fockobj.bket('X')
  %
  % There are all kinds of operations defined on Fock space objects, most 
  % importantly transposition
  %
  %   Y = X'
  %
  % scalar multiplication and direct sum
  %
  %   Y = a*X1 + b*X2
  %
  % matrix product
  %
  %   X = X1*X2
  %
  % and tensor product
  %
  %   Y = [X1 X2]
  %
  % Check out the 'example_*.mlx' live scrips to get started.
  %
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %   Peter beim Graben, BTU Cottbus-Senftenberg
  %
  % See also fockbasis, fock
  
  % == Toolbox-accessible constants ==
  properties(Constant,Access={?fock,?fockbasis})
    ROW = 1; % Column index of row basis vector identifiers in data property.
    COL = 2; % Column index of column basis vector identifiers in data property.
    VAL = 3; % Column index of values in data property
  end
  
  % == Public properties ==
  properties
    
    % Human-readable descriptive name; only used for display.
    name
    
  end

  % == Toolbox-accessible properties==
  properties(SetAccess={?fock,?fockbasis,?fockpst,?fockfvr})
  
    % Cell array containing the sparse matrix representing this object.
    data

  end
  
  % == Instance getters ==
  methods(Static)

    function obj=bket(id)
      % Creates a fockobj (unit ket) from a basis vector identifier.
      %
      %   obj=fockobj.bket(bvecid)
      %
      % arguments:
      %   id - A basis vector identifier. Use '§' as tensor product separator.
      %        For example,
      %              
      %          fockobj.bket('X§Y').
      %
      %        equals
      %
      %          kron(fockobj.bket('X'),fockobj.bket('Y'))
      %
      % returns:
      %   obj - The basis vector as a ket.
      if nargin==0                                                              % No arguments >>
        id = fock.C_BVAC;                                                       %   Use vacuum basis vector id.
      end                                                                       % <<
      id = char(id);                                                            % Convert basis vector id. to char
      id = regexprep(id,fock.C_TDELA,fock.C_TDEL);                              % Replace ASCII tensor prod. delimiter
      obj = fockobj();                                                          % Create fockobj
      obj.data = { id fock.C_BVAC 1 };                                          % Write data
      obj = obj.opretConvert();                                                 % Normalize scalar return values
    end
    
  end
  
  % == Public methods ==
  methods
    
    function t=getType(obj)
      % Determines the type of a fockobj.
      %
      %   t=obj.getType()
      %
      % arguments:
      %   obj - The fockobj.
      %   
      % returns:
      %   t   - The type, one of the following:
      %         * fock.OBJ_NULL   if the object is a null (empty),
      %         * fock.OBJ_SCALAR if the object is a scalar,
      %         * fock.OBJ_KET    if the object is a ket,
      %         * fock.OBJ_BRA    if the object is a bra, or
      %         * fock.OBJ_LOP    if the object is a linear operator.
      t = fock.OBJ_INVALID;                                                     % Default result is OBJ_INVALID
      if size(obj.data,1)==0                                                    % No data >>
        t = fock.OBJ_NULL;                                                      %   Null object (should not happen)
        return;                                                                 %   I have ready...
      end                                                                       % <<
      rbids = obj.data(:,fockobj.ROW);                                          % Get row basis vector identifiers
      cbids = obj.data(:,fockobj.COL);                                          % Get column basis vector identifiers
      if all(strcmp(rbids,string(fock.C_BVAC)))                                 % All rows vacuum >>
        if all(strcmp(cbids,string(fock.C_BVAC)))                               %   All columns vacuum, too >>
          if size(obj.data,1)==1                                                %     One data element >>
            t = fock.OBJ_SCALAR;                                                %       This is a scalar
          end                                                                   %     <<
        elseif ~any(strcmp(cbids,string(fock.C_BVAC)))                          %   << No columns vacuum >>
          t = fock.OBJ_BRA;                                                     %     This is a bra
        end                                                                     %   <<
      elseif ~any(strcmp(rbids,string(fock.C_BVAC)))                            % << No rows vacuum >>
        if all(strcmp(cbids,string(fock.C_BVAC)))                               %   All columns vacuum >>
          t = fock.OBJ_KET;                                                     %     This is a ket
        elseif ~any(strcmp(cbids,string(fock.C_BVAC)))                          %   << No columns vacuum >>
          t = fock.OBJ_LOP;                                                     %     This is linear operator
        end                                                                     %   <<
      end                                                                       % << (Everything else is flimflam)
    end

  end
  
  % == Operators ==
  % https://de.mathworks.com/help/matlab/matlab_oop/implementing-operators-for-your-class.html
  methods

    function r=eq(a,b)
      % Equality (overloads eq operator).
      %
      %    r=eq(a,b)
      %    r=(a==b)
      %
      % arguments:
      %   a - First operator, a Fock space object.
      %   b - Second operator, a Fock space object.
      %
      % returns:
      %   r - True, if a and b are equal, false otherwise

      if ~isa(a,'fockobj')||~isa(b,'fockobj')                                   % Not both fockobjs >>
        r = 0;                                                                  %   Not equal
      else                                                                      % <<
        [~,a,b]=fockobj.prepareMatrixOp(a,b);                                   %   Get matrices for fockobjs
        r = isequal(a,b);                                                       %   Compare matrices
      end                                                                       % <<
    end

    function r=ne(a,b)
      % Inequality (overloads ne operator).
      %
      %    r=ne(a,b)
      %    r=(a~=b)
      %
      % arguments:
      %   a - First operator, a Fock space object.
      %   b - Second operator, a Fock space object.
      %
      % returns:
      %   r - True, if a and b are not equal, false otherwise

      r = ~eq(a,b);                                                             % Test inequality
    end

    function r=ctranspose(a)
      % Complex conjugate transpose (overloads ' operator).
      %
      %   r=transpose(a)
      %   r=a'
      %
      % arguments:
      %   a - A Fock space object or an array of Fock space objects.
      %
      % returns:
      %   The transposed object or, if the argument is a array, the
      %   transposed array.
      %
      % throws exception:
      %   - If the argument is a non-numeric scalar but not of type
      %     'fockobj'.
      if any(size(a)>1)                                                         % Matrix of fockobj >>
        r = transpose(a);                                                       %   Transpose object matrix
        return;                                                                 %   And that was it...
      end                                                                       % <<
      r = a;                                                                    % Result gets input
      if (isnumeric(a)); return; end                                            % Transpose of numeric scalar -> ready
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      a.check();                                                                % Check a      
      tmp = r.data(:,fockobj.COL);                                              % Temp. array <- column b. vec. ids.
      r.data(:,fockobj.COL) = r.data(:,fockobj.ROW);                            % Column b.vec.ids. <- row b.vec.ids.
      r.data(:,fockobj.ROW) = tmp;                                              % Row b. vec. ids. <- temp. array
      for i=1:size(r.data,1)                                                    % Loop over cells >>
        r.data{i,fockobj.VAL} = conj(r.data{i,fockobj.VAL});                    %   Complex conjugation of values
      end                                                                       % <<
    end
    
    function r=uminus(a)
      % Unary minus (overloads unary - operator).
      %
      %   r=uminus(a)
      %   r=-a
      %
      % arguments:
      %   a - A Fock space object.
      %
      % returns:
      %   A Fock space object with all values negated.
      %
      % throws exception:
      %   - If the argument is a non-numeric scalar but not of type
      %     'fockobj'.
      if any(size(a)>1)                                                         % Matrix of fockobj >>
        r = transpose(a);                                                       %   Transpose object matrix
        return;                                                                 %   And that was it...
      end                                                                       % <<
      r = a;                                                                    % Result gets input
      if (isnumeric(a)); return; end                                            % Transpose of numeric scalar -> ready
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      a.check();                                                                % Check a      
      for i=1:size(r.data,1)                                                    % Loop over cells >>
        r.data{i,fockobj.VAL} = -r.data{i,fockobj.VAL};                         %   Negate values
      end                                                                       % <<
    end
    
    function r=plus(a,b)
      % Direct sum (overloads '+' operator).
      %
      %   r=plus(a,b)
      %   r=a+b
      %
      % arguments:
      %   a - First operator, a Fock space object or 0.
      %   b - Second operator, a Fock space object or 0.
      %
      % returns:
      %   r - The direct (aka. orthogonal) sum of the operators.
      %
      % throws exception:
      %   - If either operator is not 0 or a scalar of type 'fockobj'.
      %   - If a and b are Fock space objects and their dimensions are not 
      %     identical.

      % Initialize and check                                                    % -------------------------------------
      if isnumeric(a)&&isscalar(a)&&(a==0); a = fockobj(); end                  % a is numeric zero -> null object
      if isnumeric(b)&&isscalar(b)&&(b==0); b = fockobj(); end                  % b is numeric zero -> null object
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      assert(isa(b,'fockobj'), ...                                              % b must be a fockobj
        sprintf(fock.ERR_BADARG,'b','Must be of type ''fockobj''.'));           % ...
      if a.getType()==fock.OBJ_NULL; r = b.opretConvert(); return; end          % a is null object -> return b
      if b.getType()==fock.OBJ_NULL; r = a.opretConvert(); return; end          % b is null object -> return a

      % Prepare for matrix operation                                            % -------------------------------------
      [H,a,b]=fockobj.prepareMatrixOp(a,b);                                     % Get matrices for fockobjs
      assert(all(size(a)==size(b)), ...                                         % Matrix dimensions must be identical
        sprintf(fock.ERR_BADARG,'a or b','Must have equal dimensions.'));       % ...
      
      % Do matrix operation                                                     % -------------------------------------
      r = a+b;                                                                  % Add matrices

      % Finish result                                                           % -------------------------------------
      r = fockobj.finishMatrixOp(r,H);                                          % Get fockobj for result matrix
      r = r.opretConvert();                                                     % Normalize scalar return values
    end
    
    function r=minus(a,b)
      % Direct sum a+(-b) (overloads '-' operator).
      %
      %   r=minus(a,b)
      %   r=a-b
      %
      % arguments:
      %   a - First operator, a Fock space object or 0.
      %   b - Second operator, a Fock space object or 0.
      %
      % returns:
      %   r - The direct (aka. orthogonal) sum a+(-b).
      %
      % throws exception:
      %   - If either operator is not 0 or a scalar of type 'fockobj'.
      %   - If a and b are Fock space objects and their dimensions are not 
      %     identical.
      r = a+(-b);
    end
    
    function r=mtimes(a,b)
      % Matrix product (overloads '*' operator).
      % 
      %   r=mtimes(a,b)
      %   r=a*b
      %
      % arguments:
      %   a - First operator, a Fock space object or a numeric scalar.
      %   b - Second operator, a Fock space object or a numeric scalar.
      %
      % returns:
      %   r - The matrix product of the operators.
      %
      % throws exception:
      %   - If either operator is not a numric scalar or a scalar of type 
      %     'fockobj'.
      %   - If a and b are Fock space objects and their dimensions do not match. 

      % Initialize and check                                                    % -------------------------------------
      if isnumeric(a)&&isscalar(a); a = fockobj(a); end                         % a is numeric scalar -> scalar fockobj
      if isnumeric(b)&&isscalar(b); b = fockobj(b); end                         % b is numeric scalar -> scalar fockobj
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      assert(isa(b,'fockobj'), ...                                              % b must be a fockobj
        sprintf(fock.ERR_BADARG,'b','Must be of type ''fockobj''.'));           % ...
      if a.getType()==fock.OBJ_NULL; r = a.opretConvert(); return; end          % a is null object -> return 0
      if b.getType()==fock.OBJ_NULL; r = b.opretConvert(); return; end          % a is null object -> return 0

      % Prepare for matrix operation                                            % -------------------------------------
      [H,a,b]=fockobj.prepareMatrixOp(a,b);                                     % Get matrices for fockobjs
      if size(a,1)==1 && size(a,2)>1                                            % a is a bra >>
        if size(b,1)==1 && size(b,2)>1                                          %   b is a bra >>
          error(fock.ERR,fock.ERR_BADMULT,'<bra|*<bra|','Did yo mean a.*b?');   %     Error
        end                                                                     %   <<
      elseif size(a,1)>1 && size(a,2)==1                                        % << a is a ket >>
        if size(b,1)>1 && size(b,2)==1                                          %   b is a ket >>
          error(fock.ERR,fock.ERR_BADMULT,'|ket>*|ket>','Did yo mean a.*b?');   %     Error
        elseif size(b,1)>1 && size(b,2)>1                                       %   << b is an operator >>
          error(fock.ERR,fock.ERR_BADMULT,'|ket>*OP','');                       %     Error
        end                                                                     %   <<
      elseif size(a,1)>1 && size(a,2)>1                                         % << a is an operator >>
        if size(b,1)==1 && size(b,2)>1                                          %   b is a bra >>
          error(fock.ERR,fock.ERR_BADMULT,'OP*<bra|','');                       %     Error
        end                                                                     %  <<
      end                                                                       % <<

      % Do matrix operation                                                     % -------------------------------------
      r = a*b;                                                                  % Multiply matrices

      % Finish result                                                           % -------------------------------------
      r = fockobj.finishMatrixOp(r,H);                                          % Get fockobj for result matrix
      r = r.opretConvert();                                                     % Normalize scalar return values
    end

    function r=times(a,b)
      % Element-wise product (overloads '.*' operator).
      %
      %   r=times(a,b)
      %   r=a.*b
      %
      % arguments:
      %   a - First operator, a Fock space object or a numeric scalar.
      %   b - Second operator, a Fock space object or a numeric scalar.
      %
      % returns:
      %   r - The element-wise product of the operators.
      %
      % throws exception:
      %   - If either operator is not a numric scalar or a scalar of type 
      %     'fockobj'.
      %   - If a and b are Fock space objects and their dimensions are not 
      %     identical.
      if isnumeric(a)&&isscalar(a)||isnumeric(b)&&isscalar(b)                   % Either arg. is a numeric scalar >>
        r = mtimes(a,b);                                                        %   Invoke matrix product
      else                                                                      % << Both args. non-numeric >>
        
        % Initialize and check                                                  %   -----------------------------------
        assert(isa(a,'fockobj'), ...                                            %   a must be a fockobj
          sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));         %   ...
        assert(isa(b,'fockobj'), ...                                            %   b must be a fockobj
          sprintf(fock.ERR_BADARG,'b','Must be of type ''fockobj''.'));         %   ...

        % Prepare for matrix operation                                          %   -----------------------------------
        [H,a,b]=fockobj.prepareMatrixOp(a,b);                                   %   Get matrices for fockobjs
        assert(all(size(a)==size(b)), ...                                       %   Matrix dimensions must be identical
          sprintf(fock.ERR_BADARG,'a or b','Must have equal dimensions.'));     %   ...

        % Do matrix operation                                                   %   -----------------------------------
        r = a.*b;                                                               %   Multiply matrices element-wise

        % Finish result                                                         %   -----------------------------------
        r = fockobj.finishMatrixOp(r,H);                                        %   Get fockobj for result matrix
        r = r.opretConvert();                                                   %   Normalize scalar return values
      
      end                                                                       % <<
    end

    function r=mpower(a,b)
      % Matrix power (overloads '^' operator).
      %
      %   r=mpower(a,b)
      %   r=a^b
      %
      % arguments:
      %   a - A linear operator.
      %   b - An integer.
      %
      % returns:
      %   r - The matrix power.
      %
      % throws exception:
      %   - If a is not a linear operator or b is not an integer.
      % 
      % See also mpowerf
      
      assert(isscalar(a)&&isa(a,'fockobj'),fock.ERR_BADARG, ...                 % Arg. a must be one fockobj
        'a', [ 'Must be a scalar of type ''fockobj''. Did you mean ' ...        % ... 
               fock.HREF('fockobj/mpowerf') '?' ]);                             % ...
      assert(isscalar(b)&&~mod(b,1),fock.ERR_BADARG, ...                        % Arg. b must be one integer
        'b', [ 'Must be an integer scalar. Did you mean ' ...                   % ... 
               fock.HREF('fockobj/mpowerf') '?' ]);                             % ...
      r = mpowerf(a,b);                                                         % Invoke general matrix power
    end
    
    function r=horzcat(a,varargin)
      % Kronecker tensor product (overloads '[...]' operator).
      %
      %   r=horzcat(a1,...,an)
      %   r=[a1 ... an]
      %
      % arguments:
      %   a1      - First operator.
      %   ..., an - Further operators.
      %
      % returns:
      %   r - The Kronecker tensor product of the operators.
      %
      % throws exception:
      %   - If either operator is not a scalar of type 'fockobj'.
      assert(isscalar(a) && (isa(a,'fockobj') || isnumeric(a)), ...             % Argument #1 is bad
        '#1','Must be a numeric or ''fockobj'' scalar.');                       % ...
      r = a;                                                                    % Result gets a
      for i=1:length(varargin)                                                  % Loop over remaining arguments >> 
        b = varargin{i};                                                        %   Get current argument into b
        assert(isscalar(b) && (isa(b,'fockobj') || isnumeric(b)), ...           % Argument #i is bad
          sprintf('#%d',i),'Must be a numeric or ''fockobj'' scalar.');         % ...
        r = kron(r,b);                                                          %   Cascade invocations of kron
      end                                                                       % <<
    end

    function r=inv(a)
      % Inverts a linear operator.
      %
      %   r = inv(a)
      %
      % arguments:
      %   a - The operator.
      %
      % returns:
      %   r - The inverted operator such that a*r=eye.
      %
      % throws exception:
      %   - If a is not a scalar or linear operator.

      % Prepare for matrix operation                                            % -------------------------------------
      [H,a]=fockobj.prepareMatrixOp(a);                                         % Get matrix for fockobj
      if size(a,1)~=size(a,2)                                                   % Not a scalar or linear operator >>
        error(fock.ERR,fock.ERR_BADINV);                                        %   Error
      end                                                                       % <<
      N = size(a,1);                                                            % Square dimension
      vacrem = N && a(1,1)==0;                                                  % Vacuum basis vector to be removed.
      if vacrem; a=a(2:N,2:N); end                                              % Remove vacuum row and column

      % Do matrix operation                                                     % -------------------------------------
      r = inv(a);                                                               % Invert matrix

      % Finish result                                                           % -------------------------------------
      if vacrem                                                                 % Vacuum basis vector was removed >>
        r = [ zeros(1,N-1); r];                                                 %   Insert heading row of zeros into r
        r = [ zeros(N,1) r];                                                    %   Insert heading col.of zeros into r
      end                                                                       % <<
      r = fockobj.finishMatrixOp(r,H);                                          % Get fockobj for result matrix
      r = r.opretConvert();                                                     % Normalize scalar return values
    end
    
    function r=mpowerf(a,b)
      % Matrix power supporting matrix exponentials and non-integer powers.
      %
      %   r=mpowerf(a,b)
      %
      % arguments:
      %   a - First operator.
      %   b - Second operator.
      %
      %   One of the operators must be a numeric scalar and the other one a
      %   linear operator.
      %   
      %   NOTE: If a is numeric or b is not an integer, the function realizes 
      %         the fockobj--b or a, respectively--as a full(!) matrix.
      %
      % returns:
      %   r - The matrix power.
      %
      % throws exception:
      %   - If the operators are invalid (see above).

      % Initialize and check                                                    % -------------------------------------
      assert(isscalar(a),fock.ERR_BADARG,'a','Must not be an array');           % Check that a is not an array
      assert(isscalar(a),fock.ERR_BADARG,'b','Must not be an array');           % Check that a is not an array
      if isa(a,'fockobj')                                                       % a is a fockobj >>
        a = a.opretConvert();                                                   %   See if it is actually a scalar
      else                                                                      % << a is not a fockobj >>
        assert(isnumeric(a),fock.ERR_BADARG,'b', ...                            %   Then it must be numeric
          'Must be numeric or of type ''fockobj''.');                           %   ...
      end                                                                       % <<
      if isa(b,'fockobj')                                                       % b is a fockobj >>
        b = b.opretConvert();                                                   %  See if it is actually a scalar
      else                                                                      % << b is not a fockobj >>
        assert(isnumeric(b),fock.ERR_BADARG,'b', ...                            %   Then it must be numeric
          'Must be numeric or of type ''fockobj''.');                           %   ...
      end                                                                       % <<
      assert(~isa(a,'fockobj')||~isa(b,'fockobj'), ...                          % At least argument must be numeric
        fock.ERR_BADARG,'a or b','At least one must be numeric.');              % ...

      % Prepare for matrix operation                                            % -------------------------------------
      if isnumeric(a)&&isnumeric(b); r = a^b; return; end                       % Two scalars -> just do it and return
      if isa(a,'fockobj')                                                       % The (only!) fockobj is a >>
        [H,a]=fockobj.prepareMatrixOp(a);                                       %   Get matrix for a
        assert(size(a,1)==size(a,2),fock.ERR_BADARG, ...                        %   Check if a is a linear operator
          'a','Matrix power only defined for scalars and linear operators');    %   ...
        N = size(a,1);                                                          %   Get square dimension
        a = a(2:N,2:N);                                                         %   Remove vacuum row and column
        if mod(b,1)                                                             %   b not integer >> 
          if fock.CHECK; warning('Expaning argument a to full matrix!'); end    %     Warning in CHECK mode
          a=full(a);                                                            %     Convert a to full matrix
        end                                                                     %   <<
      else                                                                      % << The (only!) fockobj is b >>
        [H,b]=fockobj.prepareMatrixOp(b);                                       %   Get matrix for b
        assert(size(b,1)==size(b,2),fock.ERR_BADARG, ...                        %   Check if a is a linear operator
          'b','Matrix power only defined for scalars and linear operators');    %   ...
        N = size(b,1);                                                          %   Get square dimension
        b = b(2:N,2:N);                                                         %   Remove vacuum row and column
        if fock.CHECK; warning('Expaning argument b to full matrix!'); end      %   Warning in CHECK mode
        b = full(b);                                                            %   Exponential --> full matrix b!
      end                                                                       %  <<

      % Do matrix operation                                                     % -------------------------------------
      r = a^b;                                                                  % Matrix power

      % Finish result                                                           % -------------------------------------
      r = [ zeros(1,N-1); r];                                                   % Restore vacuum row
      r = [ zeros(N,1) r];                                                      % Restore vacuum column
      r = fockobj.finishMatrixOp(r,H);                                          % Get fockobj for result matrix
      r = r.opretConvert();                                                     % Normalize scalar return values

    end
    
    function r=kron(a,b)
      % Kronecker tensor product.
      %
      %   r=kron(a,b)
      %
      % arguments:
      %   a - First operator.
      %   b - Second operator.
      %
      % returns:
      %   r - The Kronecker tensor product of the operators.
      %
      % throws exception:
      %   - If either operator is not a scalar of type 'fockobj'.
      if isnumeric(a)&&isscalar(a); a = fockobj(a); end                         % Convert numeric -> fockobj sclar
      if isnumeric(b)&&isscalar(b); b = fockobj(b); end                         % Convert numeric -> fockobj sclar
      assert(isa(a,'fockobj'),fock.ERR_BADARG,'a', ...                          % a is bad
        'Must be a numeric or ''fockobj'' scalar.');                            % ...
      assert(isa(b,'fockobj'),fock.ERR_BADARG,'b', ...                          % b is bad
        'Must be a numeric or ''fockobj'' scalar.');                            % ...
      na = size(a.data,1);                                                      % Get number of elements in a
      nb = size(b.data,1);                                                      % Get number of elements in b
      r = fockobj();                                                            % Create result fockobj
      r.data = cell(na*nb,3);                                                   % Pre-allocate data
      for i=1:na                                                                % Loop over elements in a >>
        ra = a.data{i,fockobj.ROW};                                             %   Get row basis vector id. in a
        ca = a.data{i,fockobj.COL};                                             %   Get column basis vector id. in a
        va = a.data{i,fockobj.VAL};                                             %   Get value in a
        for j=1:nb                                                              %   Loop over elements in b >>
          rb = b.data{j,fockobj.ROW};                                           %     Get row basis vector id. in b
          cb = b.data{j,fockobj.COL};                                           %     Get col. basis vector id. in b
          vb = b.data{j,fockobj.VAL};                                           %     Get value in b
          k  = (i-1)*nb+j;                                                      %     Compute element index in result
          r.data{k,fockobj.ROW} = [ ra fock.C_TDEL rb ];                        %     Store combined row b. vec. id.
          r.data{k,fockobj.COL} = [ ca fock.C_TDEL cb ];                        %     Store combined column b. vec. id.
          r.data{k,fockobj.VAL} = va*vb;                                        %     Store product of values
        end                                                                     %   <<
      end                                                                       % <<
      r.data(:,1) = strrep(r.data(:,1),[ fock.C_TDEL fock.C_BVAC ],'');         % Remove \otimes\emptyset from row...
      r.data(:,2) = strrep(r.data(:,2),[ fock.C_TDEL fock.C_BVAC ],'');         % ... and column b. vec. ids. of result
      r.data(:,1) = strrep(r.data(:,1),[ fock.C_BVAC fock.C_TDEL ],'');         % Remove \emptyset\otimes from row...
      r.data(:,2) = strrep(r.data(:,2),[ fock.C_BVAC fock.C_TDEL ],'');         % ... and column b. vec. ids. of result
      r = r.opretConvert();                                                     % Normalize scalar return values
    end

    function r=real(a)
      % Real part of complex values.
      %
      %   r=real(a)
      %
      % arguments:
      %   a - A Fock space object.
      %
      % returns:
      %   A fock space object whose values are the real parts of the
      %   argument's values.
      %
      % throws exception:
      %   - If the argument is neither a numeric scalar nor of type 'fockobj'.
      a.check();                                                                % Check argument
      r = a;                                                                    % Copy argument to result
      if (isnumeric(a)); r = real(a); return; end                               % Numeric scalar -> return real part
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      for i=1:size(r.data,1)                                                    % Loop over cells >>
        r.data{i,fockobj.VAL} = real(r.data{i,fockobj.VAL});                    %   Real part of values
      end                                                                       % <<
      r = r.removeZeros();                                                      % Remove zero-entries
      r = r.opretConvert();                                                     % Adjust object type
    end

    function r=imag(a)
      % Imaginary part of complex values.
      %
      %   r=imag(a)
      %
      % arguments:
      %   a - A Fock space object.
      %
      % returns:
      %   A fock space object whose values are the imaginary parts of the
      %   argument's values.
      %
      % throws exception:
      %   - If the argument is neither a numeric scalar nor of type 'fockobj'.
      a.check();                                                                % Check argument
      r = a;                                                                    % Copy argument to result
      if (isnumeric(a)); r = imag(a); return; end                               % Numeric scalar -> return imag. part
      assert(isa(a,'fockobj'), ...                                              % a must be a fockobj
        sprintf(fock.ERR_BADARG,'a','Must be of type ''fockobj''.'));           % ...
      for i=1:size(r.data,1)                                                    % Loop over cells >>
        r.data{i,fockobj.VAL} = imag(r.data{i,fockobj.VAL});                    %   Imaginary part of values
      end                                                                       % <<
      r = r.removeZeros();                                                      % Remove zero-entries
      r = r.opretConvert();                                                     % Adjust object type
    end

  end

  % == Protected methods ==
  methods(Access=protected)

    function displayScalarObject(obj)
      % Custom display method (overload).
      obj.check();                                                              % Check object
      s = string(obj.name);                                                     % Get name property as string
      if size(obj.data,1)==0                                                    % Null object >>
        fprintf('  null %s\n',fock.HREF('fockobj'));                            %   Display general info
      else                                                                      % << Non-null object >>
        fprintf('  %s\n',fock.HREF('fockobj'));                                 %   Display general info
      end                                                                       % <<      
      if ~strcmp(s,"")                                                          % Non-empty name
        fprintf('\n    Name: <strong>%s</strong> ',obj.name);                   %   Display name
      end                                                                       % <<
      if size(obj.data,1)==0; return; end                                       % Null object -> all done
      fprintf('\n    Type: ');                                                  % Display object type
      switch obj.getType()                                                      % ...
        case fock.OBJ_SCALAR; fprintf('<strong>scalar</strong> ');              % ...
        case fock.OBJ_KET;    fprintf('<strong>ket</strong> ');                 % ...
        case fock.OBJ_BRA;    fprintf('<strong>bra</strong> ');                 % ...
        case fock.OBJ_LOP;    fprintf('<strong>linear operator</strong> ');     % ...
        otherwise;            fprintf(2,'<strong>invalid </strong> ');          % ... (should not happen!)
      end                                                                       % ...
      fprintf('\n\n');                                                          % ...
      for i=1:size(obj.data,1)                                                  % Loop over data elements >>
        s1 = obj.data{i,fockobj.ROW};                                           %   Get row basis vector identifier
        s2 = obj.data{i,fockobj.COL};                                           %   Get column basis vector identifier
        v  = obj.data{i,fockobj.VAL};                                           %   Get value
        s  = [ '(' s1 ',' s2 ')' ];                                             %   Make element identifier
        fprintf('    %-20s\t%f',s,real(v));                                     %   Display real part of value
        if ~isreal(v)                                                           %   Value has imaginary part >>
          if imag(v)<0; fprintf(' - '); else; fprintf(' + '); end               %     Display sign of imag. part
          fprintf('%fi',abs(imag(v)));                                          %     Display abs. imag. value
        end                                                                     %   <<
        fprintf('\n');                                                          %   ...
      end                                                                       % <<
    end
    
  end
  
  % == Toolbox-accessible constructors ==
  methods(Access={?fock,?fockbasis,?fockpst,?fockfvr})

    function check(obj)
      % In-depth consistency check.
      % The method does nothing if fock.CHECK is false.
      %
      %    check(obj)
      %
      % arguments:
      %    obj - The Fock space object to check.
      %
      % throws exception:
      % - If the check fails.
      if ~fock.CHECK; return; end                                               % Return if checking is disabled
      if isempty(obj.data); return; end                                         % No data -> ok
      assert(iscell(obj.data)&&size(obj.data,2)==3,...                          % Check type and dimesions of obj.data
        fock.ERR_CHECK,'obj.data must be a Nx3 cell array');                    % ...
      cids = strcat(obj.data(:,1),fock.C_RCDEL,obj.data(:,2));                  % Get matrix cell identifiers
      assert(length(unique(cids))==size(obj.data,1), ...                        % Check if cell idenfiers are unique
        fock.ERR_CHECK,'obj.data containes duplicate entries');                 % ...
      if any(cell2mat(obj.data(:,3))==0)                                        % Look for zero values, if any >>
        if ~obj.getType()==fock.OBJ_SCALAR                                      %   If not a scalar >>
          warning(fock.ERR_CHECK,'obj.data contains zero values');              %     Warning
        end                                                                     %   <<
      end                                                                       % <<
    end

    function obj=fockobj(arg1,arg2,arg3)
      % Constructor create Fock space object.
      %
      %   obj=fockobj()
      %   obj=fockobj(x)
      %   obj=fockobj(R,C,V)
      %
      % arguments:
      %   x - A numeric scalar.
      %   R - A cell column vector containing row basis vector identifiers.
      %   C - A cell column vector containing column basis vector identifiers.
      %   V = A cell column vector containing numeric scalars.
      %
      % returns:
      %   The newly created Fock space object.
      %
      % throws exception:
      % - If invoked with two or more than three arguments.
      % - If n is not a numeric scalar.
      % - If R, C, or V are not of the format specified above.
      assert(~verLessThan('matlab','9.4'),fock.ERR_MVERSION);                   % Check Matlab version
      switch nargin                                                             % Branch for number of arguments >>
        case 0                                                                  %   0: Create null object
          % Nothing to be done                                                  %  
        case 1                                                                  %   1: Create scalar object
          assert(isscalar(arg1)&&isnumeric(arg1), ...                           %     arg1 must be a numeric scalar
            fock.ERR_BADARG,'arg1','Must be a numeric scalar.');                %     ...
          obj.data = { fock.C_BVAC fock.C_BVAC arg1 };                          %     Store data
        case 3                                                                  %   3: Create from raw data
          obj.data = [arg1 arg2 arg3];                                          %     Store date
        otherwise                                                               %   All other numbers of arguments:
          error(fock.ERR_BADNUMARG,'Needs 0, 1, or 3 arguments.');              %     Not allowed
      end                                                                       % <<
      obj.check();                                                              % Check newly created object.
    end

  end
  
  % == Private methods ==
  methods(Access=private)
  
    function r=opretConvert(obj)
      % Converts a fockobj to be returned from an operator.
      % The method will return
      % - 0 if obj is of type fock.OBJ_NULL.
      % - a numeric scalar if obj is of type fock.OBJ_SCALAR.
      r = obj;                                                                  % Default return value is the object
      n = size(obj.data,1);                                                     % Get number of data elements
      if n==0                                                                   % No elements >>
        r = 0;                                                                  %   Return scalar zero
      elseif n==1                                                               % << One element >>
        if all(obj.data{1,fockobj.ROW}==fock.C_BVAC) ...                        %   Object is a scalar >>
            && all(obj.data{1,fockobj.COL}==fock.C_BVAC)                        %   ...
          r = obj.data{1,fockobj.VAL};                                          %     Return scalar element value
        end                                                                     %   <<
      end                                                                       % <<
    end
    
    function r=removeZeros(obj)
      % Removes data entries with zero values.
      obj.check();                                                              % Check argument
      r = obj;                                                                  % Copy argument to return value
      f = [r.data{:,fockobj.VAL}]==0;                                           % Find zero values
      r.data(f,:)=[];                                                           % Remove found rows from data
    end
    
  end
  
  % == Private static methods ==
  methods(Static,Access=private)
    
    function [H,a,b]=prepareMatrixOp(a,b)
      % Creates a Fock subspace basis enclosing the argument(s) and realizes 
      % them.
      % 
      %   [H,a]=prepareMatrixOp(a)
      %   [H,a,b]=prepareMatrixOp(a,b)
      %
      % arguments:
      %   a - First Fock space object.
      %   b - Second Fock space object, optional.
      %
      % returns:
      %   a - A sparse matrix realizing input a in Fock space H.
      %   b - A sparse matrix realizing input b in Fock space H.
      %   H - The Fock subspace basis enclosing a and b.
      %
      % See also finishMatrixOp
      V = cell(nargin,1);                                                       % Pre-allocate cell vector of args.
      a.check();                                                                % Check first argument
      V{1} = a;                                                                 % Store first argument in V
      if nargin>1                                                               % Two arguments >>
        b.check();                                                              %   Check second argument
        V{2} = b;                                                               %   Store second argument in V
      end                                                                       % <<
      H = fockbasis(V);                                                         % Get Fock space enclosing all args.
      a = H.realize(a);                                                         % Realize a in enclosing space
      if nargin>1                                                               % Two arguments >>
        b = H.realize(b);                                                       %   Realize b in enclosing space
      end                                                                       % <<
    end
    
    function r=finishMatrixOp(r,H)
      % Unrealizes a sparse matrix from a Fock space.
      %
      %   r=finishMatrixOp(r,H)
      %
      % arguments:
      %   r - A sparse matrix representing a Fock space object in H.
      %   H - The Fock subspace basis in which input r is realized.
      %
      % returns:
      %   r - A Fock space object unrealizing input r from Fock space H.
      %
      % See also prepareMatrixOp
      r = H.unrealize(r);                                                       % Unrealize r from H
      if isa(r,'fockobj'); r.check(); end                                       % Check result
    end
  
  end
  
end

% EOF