classdef fockbasis < matlab.mixin.CustomDisplay
  % A Fock basis defines a set of orthonormal basis vectors, each bearing a 
  % unique identifier. Fock bases can be added (vector sum)
  % 
  %   F = F1 + F2, 
  %
  % where the vector sum is identical with the orthogonal sum iff the bases
  % of F1 and F2 are disjoint, and multiplied (tensor product) 
  %
  %   F = [F1 F2]
  %
  % Check out the 'example_*.mlx' live scrips to get started.
  %
  % authors:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %   Peter beim Graben, BTU Cottbus-Senftenberg
  %
  % See also fockobj, fock
  
  % == Write-protected properties ==
  properties(Access=private)

    % Cell vector of basis vector identifiers.
    bid
    
    % Map of basis vector identifiers to basis indexes.
    bidm

    % Map of basis indexes to basis vector identifiers.
    bvec
    
  end
  
  % == Public constructors and instance getters ==
  methods
    
    function obj=fockbasis(V)
      % Creates a Fock subspace basis.
      %
      %   obj=fockbasis(V)
      %
      % arguments:
      %   V   - A vector of basis vector identifiers or a vector of fockobj.
      %
      % returns:
      %   obj - The Fock space.
      %
      % throws exception: 
      %   - If the argument is not vector.
      assert(~verLessThan('matlab','9.4'),fock.ERR_MVERSION);                   % Check Matlab version
      obj.bidm = containers.Map();                                              % Create basis vector ids. -> index map
      obj.bvec = containers.Map();                                              % Create basis vector map
      switch nargin                                                             % Switch for number of input arguments
        case 0                                                                  %   None:
          obj.bidm = containers.Map(fock.C_BVAC,0);                             %     Create vacuum unit state
        case 1                                                                  %   One:
          assert(isvector(V),fock.ERR_BADARG,'A','Must be a vector.');          %     Argument must be a vector
          A = fock.vec2cell(V,'fockobj');                                       %     Convert to fockobj vector
          if ~isempty(A)                                                        %     Succeeded >>
            obj.bidm = containers.Map();                                        %       Create empty b.vec ids. map
            for i=1:length(A)                                                   %       Loop over fockobjs >>
              if ~isempty(A{i}.data)                                            %         If not empty >>
                B = [ A{i}.data(:,fockobj.ROW); A{i}.data(:,fockobj.COL) ];     %           Get all basis vector ids.
                obj.bidm = [ obj.bidm; containers.Map(B,zeros(length(B),1)) ];  %           Append to basis vec. ids. map
              end                                                               %         <<
            end                                                                 %       <<
          else                                                                  %     << Not a fockobj vector
            A = fock.vec2cell(V,'char');                                        %       Convert to char array vector
            if ~isempty(A)                                                      %       Succeeded >>
              A = strrep(A,fock.C_TDELA,fock.C_TDEL);                           %         Replace tensor prod.delimiter
              obj.bidm = containers.Map(A,zeros(length(A),1));                  %         Initialize base ids. map
            else                                                                %       << Not a char arr.vec.either >>
              error(fock.ERR,fock.ERR_BADARG,'A',...                            %         Error
                'Must be a vector of char or fockobj.');                        %         ...
            end                                                                 %       <<
          end                                                                   %     <<
        otherwise                                                               %   More than one argument:
          error(fock.ERR_TOOMANYARGS);                                          %     Throw error
      end                                                                       % <<
      obj = obj.reindex();                                                      % Rebuild indexes and maps
      obj.check();                                                              % Check
    end

    function r=sector(obj,varargin)
      % Retrieves sectors of the Fock basis.
      %
      %   r=obj.sector(n)
      %
      % arguments:
      %   obj - A Fock subspace basis.
      %   n   - A positive integer or a vector thereof.
      %
      % return:
      %   The subspace cataining the sectors mentioned in n. The null
      %   sector containing the vacuum basis vector is always included in the
      %   result.
      n = cell2mat(fock.vec2cell(varargin,'double'));                           % Convert input to numeric array
      assert(~isempty(n)&&isvector(n)&&all(n>0),fock.ERR_BADARG, ...            % Check if not empty and of ...
        'n','A positive integer or a vector thereof.');                         % ... positive integers
      map(n)=1;                                                                 % Create sector index map
      keys = obj.bidm.keys;                                                     % Get basis vector identifiers
      B = cell(1,length(keys));                                                 % Pre-allocate vector of b. vec. ids.
      j=1;                                                                      % Pointer into B
      for i=1:length(keys)                                                      % Loop over basis vec. ids. >>
        if string(keys{i})==fock.C_BVAC; continue; end                          %   Ignore vacuum basis vector
        l=size(strfind(keys{i},fock.C_TDEL),2)+1;                               %   Get sector index
        if l<=length(map) && map(l)                                             %   Sector to be retrieved >>
          B{j}=keys{i};                                                         %     Add basis vector id. to B
          j=j+1;                                                                %     Increment pointer
        end                                                                     %   <<
      end                                                                       % <<
      B=B(1,1:j-1);                                                             % Truncate vector of b. vec. ids.
      r=fockbasis();                                                            % Default result is a empty space
      if ~isempty(B); r=fockbasis(B); end                                       % Create sub-space
    end
 
  end
  
  % == Public methods ==
  
  methods
    
    function b=isCanonical(obj)
      % Determines if this Fock basis is canonical.
      %
      %   b=isCanonical(obj)
      %
      % arguments:
      %   obj -  A Fock subspace basis.
      %
      % returns:
      %   b   - TRUE if the basis is canonical, FALSE if it is custom.
      b = isempty(obj.bvec);                                                    % Canonical = no basis vectors 
    end
    
    function [n,nmax]=getNumSectors(obj)
      % Determines the number of Fock space sectors covered by this basis.
      %
      %   [n,name] = numSectors(obj)
      %
      % arguments:
      %   obj  -  A Fock subspace basis.
      %
      % returns:
      %   n    - The number of sectors.
      %   nmax - The greatest sector index.
      n=0;                                                                      % Max. #tensor products
      m=intmax;                                                                 % Min. #tensor products
      keys = obj.bidm.keys;                                                     % Get basis vector identifiers
      for i=1:length(keys)                                                      % Loop over basis vector ids. >>
        if string(keys{i})==fock.C_BVAC; continue; end                          %   Ignore vacuum basis vector
        l=size(strfind(keys{i},fock.C_TDEL),2);                                 %   Get #tensor products
        n=max(n,l);                                                             %   Aggregate maximum
        m=min(m,l);                                                             %   Aggregate minimum
      end                                                                       % <<
      nmax = n+1;                                                               % Store greatest sector index
      if (n==m)                                                                 % All basis vectors in one sector >>
        n=2;                                                                    %   Hilbert space + vacuum space
      else                                                                      % << Not all in one sector >>
        n=n+2;                                                                  %   #sectors = max. #tensor prods. +1
      end                                                                       % <<
    end

    function n=getDim(obj)
      % Determines the dimension of this Fock subspace basis.
      %
      %   n=getDim(obj)
      %
      % arguments:
      %   obj - A Fock subspace basis.
      %
      % returns:
      %   n   - The dimension.
      n = length(obj.bidm.keys);                                                % Return number of basis vector ids.
    end
    
    function bid=getBvecId(obj,n)
      % Retrieves a basis vector identifier.
      %
      %   bid = getBvecId(obj,n)
      %
      % arguments:
      %   obj - A Fock subspace basis.
      %   n   - The index of the basis vector, 1 <= n <= obj.getDim().
      %
      % returns:
      %   bid - The basis vector identifier.
      %
      % throws exception:
      % - If n is not an integer scalar or out of range.
      assert(isscalar(n)&&~mod(n,1));                                           % Verify that n is an integer scalar
      assert(n>=1&&n<=length(obj.bidm.keys));                                   % Verify that n is within range
      bid = obj.bid{n};                                                         % Lookup and return basis vector id.
    end
    
    function m=realize(obj,f)
      % Converts a fock space object to a sparse matrix.
      %
      %   m=obj.realize(f)
      %
      % arguments:
      %   obj - The Fock subspace basis to realize on.
      %   f   - The Fock space object to realize.
      %
      % returns:
      %   m   - A sparse matrix representing a in the basis of obj.
      %
      % Examples:
      %
      %   - Display full matrix realizing fockobj f
      %
      %     full(fockbasis(f).realize(f)) 
      %
      % See also unrealize
      assert(isa(f,'fockobj'),fock.ERR_BADARG,'f', ...                          % Argument f must be a fpckobj
        'Must be of type ''fockobj''.');                                        % ...
      obj.check();                                                              % Check
      assert(obj.isCanonical,sprintf(fock.ERR_NOTIMPL,'Custom base vectors'));  % Not implemented for custom base
      assert(all(isKey(obj.bidm,f.data(:,fockobj.ROW))),fock.ERR_BADARG,'f', ...% Check row basis vector ids.
        'Contains row base identifiers outside the space.');                    % ...
      assert(all(isKey(obj.bidm,f.data(:,fockobj.COL))),fock.ERR_BADARG,'f', ...% Check column basis vector ids.
        'Contains col base identifiers outside the space.');                    % ...
      R = cell2mat(values(obj.bidm,f.data(:,fockobj.ROW)));                     % Convert row b. vec. ids. to indexes
      C = cell2mat(values(obj.bidm,f.data(:,fockobj.COL)));                     % Convert column b. vec. ids. to indexes
      V = cell2mat(f.data(:,fockobj.VAL));                                      % Get values
      t = f.getType();                                                          % Type of Fock space object is >>
      if     t==fock.OBJ_NULL;   n=0; m=0;                                      %   OBJ_NULL       -> matrix size: (0x0)
      elseif t==fock.OBJ_SCALAR; n=1; m=1;                                      %   OBJ_SCALAR     -> matrix size: (1x1)
      elseif t==fock.OBJ_KET;    n=length(obj.bidm); m=1;                       %   OBJ_KET        -> matrix size: (Nx1)
      elseif t==fock.OBJ_BRA;    n=1; m=length(obj.bidm);                       %   OBJ_KET        -> matrix size: (1xN)
      elseif t==fock.OBJ_LOP;    n=length(obj.bidm); m=length(obj.bidm);        %   OBJ_LOP        -> matrix size: (NxN)
      else; assert(false,sprintf(fock.ERR_OBJINVAL,'a'));                       %   something else -> exception
      end                                                                       % <<
      m = sparse(R,C,V,n,m);                                                    % Create sparse matrix
    end

    function f=unrealize(obj,m)
      % Converts a (sparse) matrix to a fock space object.
      %
      %   f=fockobj.unrealize(obj,m)
      %
      % arguments:
      %   obj - The Fock subspace basis to unrealize from.
      %   m   - A sparse matrix representing f in the basis of obj.
      %
      % returns:
      %   f   - The unrealized Fock space object.
      %
      % See also realize
      obj.check();                                                              % Check
      assert(obj.isCanonical,sprintf(fock.ERR_NOTIMPL,'Custom base vectors'));  % Not implemented for custom basis
      m = sparse(m);                                                            % Convert input matrix to sparse
      [R,C,V] = find(m);                                                        % Get row/column indexes and values
      R = obj.bid(R,1);                                                         % Get row basis vector ids.
      C = obj.bid(C,1);                                                         % Get column basis vector ids.
      if isrow(V); V=V.'; end                                                   % Ensure values are a column
      V = num2cell(V);                                                          % Convert values to cell array
      f = fockobj(R,C,V);                                                       % Create Fock space object
    end
    
  end

  % == Operators ==
  % https://de.mathworks.com/help/matlab/matlab_oop/implementing-operators-for-your-class.html
  methods

    function r=plus(a,b)
      % Vector sum (overloads '+' operator). The vector sum is identical
      % with the direct (or orthogonal) sum iff the bases of the operands
      % are disjoint.
      %
      %   r=plus(a,b)
      %   r=a+b
      %
      % arguments:
      %   a - First operand, a Fock basis object.
      %   b - Second operand, a Fock basis object.
      %
      % returns:
      %   r - The sum of the operands.
      %
      % throws exception:
      %   - If either operand is not a scalar of type 'fockbasis'.
      assert(isscalar(a)&&isa(a,'fockbasis'),...                                % Check if a is a fockbasis scalar
        'a','Must be a scalar of type ''fockbasis''');                          % ...
      assert(isscalar(b)&&isa(b,'fockbasis'),...                                % Check if b is a fockbasis scalar
        'b','Must be a scalar of type ''fockbasis''');                          % ...
      a.check();                                                                % Check argument #1
      b.check();                                                                % Check argument #2
      assert(a.isCanonical()&&b.isCanonical, ...                                % Not implemented for custom basis
        fock.ERR_NOTIMPL, 'Custom basis vectors');                              % ...
      r = fockbasis([ a.bidm.keys b.bidm.keys ]);                               % Create sum space
    end

    function r=horzcat(a,varargin)
      % Kronecker tensor product (overloads '[...]' operator).
      %
      %   r=horzcat(a1,...,an)
      %   r=[a1 ... an]
      %
      % arguments:
      %   a1      - First operand.
      %   ..., an - Further operand.
      %
      % returns:
      %   r - The Kronecker tensor product of the operaand.
      %
      % throws exception:
      %   - If either operand is not a scalar of type 'fockbasis'.
      % 
      % TODO: Basis vector identifier concatenation should go into a
      %       separate function.
      assert(isscalar(a)&&isa(a,'fockbasis'), ...                               & Check if a1 is a fockbasis scalar
        '#1','Must be a scalar of type ''fockbasis''');                         % ...
      a.check();                                                                % Check argument #1
      keysa = a.bidm.keys;                                                      % Get base identifiers a
      for i=1:length(varargin)                                                  % Loop over remaining arguments >>
        b = varargin{i};                                                        %   Get next argument into b
        assert(isscalar(b)&&isa(b,'fockbasis'),sprintf('#%d',i),...             %   Check if it is a fockbasis
          'Must be a scalar of type ''fockbasis''');                            %   ...
        b.check();                                                              %   Check it
        keysb = b.bidm.keys;                                                    %   Get basis vector identifiers b
        na    = length(keysa);                                                  %   Get number of b. vec. ids. in a
        nb    = length(keysb);                                                  %   Get number of b. vec. ids. in b
        keysr = cell(1,na*nb);                                                  %   Pre-allocate resulting list
        for j=1:na                                                              %   Loop over basis vector ids. a >>
          bida = keysa{j};                                                      %     Get basis vector identifier a
          vaca = all(bida==fock.C_BVAC);                                        %     Detect vacuum basis id.
          for k=1:nb                                                            %     Loop over basis vector ids. b >>
            bidb = keysb{k};                                                    %       Get basis vector identifier b   
            vacb = all(bidb==fock.C_BVAC);                                      %       Detect vacuum basis vector id.
            if vaca && vacb                                                     %       Both vacuum >>
              bidr = fock.C_BVAC;                                               %         Result is vacuum
            elseif vaca                                                         %       << id. a is vacuum >>
              bidr = bidb;                                                      %         Result is id. b
            elseif vacb                                                         %       << id. b is vacuum >>
              bidr = bida;                                                      %         Result is id. a
            else                                                                %       << both not vacuum >>
              bidr = [ bida fock.C_TDEL bidb ];                                 %         Concatenate ids a and b.
            end                                                                 %       <<
            keysr{(j-1)*nb+k} = bidr;                                           %       Store resulting basis vec. id.
          end                                                                   %     <<
        end                                                                     %   <<
        keysa = keysr;                                                          %   Basis vector ids. a get result
      end                                                                       % <<
      r = fockbasis(keysa);                                                     % Create tensor product space
    end  

    function r=kron(a,b)
      % Kronecker tensor product.
      %
      %   r=kron(a,b)
      %
      % arguments:
      %   a - First operand.
      %   b - Second operand.
      %
      % returns:
      %   r - The Kronecker tensor product of the operands.
      %
      % throws exception:
      %   - If either operand is not a scalar of type 'fockbasis'.
      r = horzcat(a,b);                                                         % Invoke horizontal concatenation
    end

    function r=mpower(a,b)
      % Tensor product of order b (overloads '^' operator).
      %
      %   r=mpower(a,b)
      %   r=a^b
      %
      % arguments:
      %   a - First operand, scalar of type 'fockbasis'.
      %   b - Second operand, a positive integer.
      %
      % returns:
      %   r - The tensor product of order b.
      %
      % throws exception:
      %   - If first operand is not a scalar of type 'fockbasis' or second
      %     operand is not a positive integer.
      
      assert(isscalar(a)&&isa(a,'fockbasis'),fock.ERR_BADARG, ...               % Arg. a must be one fockbasis
        'a', 'Must be a scalar of type ''fockbasis''.');                        % ...
      assert(isscalar(b)&&~mod(b,1)&&b>0,fock.ERR_BADARG, ...                   % Arg. b must be one positive integer
        'b', 'Must be a positive integer scalar.');                             % ...
      r = a; for i=2:b; r = kron(r,a); end                                      % Compute b-th order tensor product
    end

  end

  % == Protected Methods ==
  methods(Access=protected)

    function check(obj)
      % In-depth consistency check.
      % The method does nothing if fock.CHECK is false.
      %
      %    check(obj)
      %
      % arguments:
      %    obj - The Fock space to check.
      %
      % throws exception:
      % - If the check fails.
      if ~fock.CHECK; return; end                                               % Return if checking is disabled
      bidm2 = containers.Map(obj.bid,1:length(obj.bid));                        % Create map of basis vector ids.
      assert(isequal(obj.bidm,bidm2),...                                        % Assert equality to internal map
        fock.ERR_CHECK,'Base identifier list and map not consistent.');         % ...
    end
    
    function obj=reindex(obj)
      % Updates the basis vector identifier to basis vector index map.
      %
      % arguments:
      %   obj - The fock space.
      %
      % returns:
      %   obj - The re-indexed Fock basis.
      if obj.bidm.isKey(fock.C_BVAC); obj.bidm.remove(fock.C_BVAC); end         % Remove vacuum basis vector
      if ~isempty(obj.bidm)                                                     % If other basis vectors exist >>
        obj.bidm = containers.Map(obj.bidm.keys,2:obj.bidm.length+1);           %   Re-index them
      end                                                                       % <<
      obj.bidm = [ containers.Map({fock.C_BVAC},{1}); obj.bidm ];               % Add vacuum basis vector w/ index #1
      
      % NOTE: This is cowardly assuming no stable order of the key set -->
      keys = obj.bidm.keys;                                                     % Fock space dimension
      obj.bid = cell(length(keys),1);                                           % Get basis vector identifiers
      for i=1:length(keys)                                                      % Loop over basis vector ids. >>
        s = keys{1,i};                                                          %   Get basis vector identifier
        j = values(obj.bidm,{s});                                               %   Get basis vector index
        obj.bid{j{1},1}=s;                                                      %   Write to basis vector id. list
      end                                                                       % <<
      % <--
    end
    
    function displayScalarObject(obj)
      % Custom display method (overload).
      obj.check();
      if isempty(obj.bidm)
        fprintf('  empty %s\n',fock.HREF('fockbasis'));
      else
        fprintf('  %s\n\n',fock.HREF('fockbasis'));
        fprintf('    Dimensions: <strong>%d</strong>\n',length(obj.bidm));
        n = obj.getNumSectors();
        fprintf('    Sectors   : <strong>%d</strong>\n',n);
        fprintf('    Type      : '); 
        if obj.isCanonical()
          fprintf('<strong>canonical</strong>\n');
        else
          fprintf('<strong>custom</strong>\n');
        end
        fprintf('\n');
        for i=1:length(obj.bid)
          fprintf('    %-20s\t%d\n',obj.bid{i},i);
        end
        %keys = obj.bidm.keys;
        %for i=1:length(keys)
        %  fprintf('    %-20s\t%d\n',keys{i},obj.bidm(keys{i}));
        %end
      end
    end
    
  end

end

% EOF