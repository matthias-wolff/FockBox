classdef qcirc
  % Library of quantum gates as Fock space operators.
  %
  % author:
  %   Matthias Wolff, BTU Cottbus-Senftenberg
  %
  % TODO:
  %   https://github.com/matthias-wolff/FockBox/issues?q=is%3Aopen+is%3Aissue+project%3Amatthias-wolff%2FFockBox%2F1
  %
  % See also fockobj
  
  %% == Constants ==
  
  properties (Constant)

    % The qubit Fock space basis { |eps>, |0>, |1> }.
    basis = fockbasis({qcirc.b0 qcirc.b1});
    
    % Identity gate operator (1 qubit).
    id = qcirc.opFromMatrix([1 0; 0 1],'Identity');

    % Pauli-X gate operator (1 qubit).
    X = qcirc.opFromMatrix([0 1; 1 0],'Pauli-X');

    % Pauli-Y gate operator (1 qubit).
    Y = qcirc.opFromMatrix([0 -1i; 1i 0],'Pauli-Y');
  
    % Pauli-Z gate operator (1 qubit).
    Z = qcirc.opFromMatrix([1 0; 0 -1],'Pauli-Z');

    % Hadamard gate operator (1 qubit).
    H = qcirc.opFromMatrix(1/sqrt(2)*[1 1; 1 -1],'Hadamard');

    % Phase gate S operator (1 qubit), S = mpowerf(qcirc.Z,1/2).
    % 
    % remarks:
    %   - The inverse is obtained by qcirc.S'
    S = qcirc.opFromMatrix([1 0; 0 1i],'S phase');

    % Phase gate T operator (1 qubit), T = mpowerf(qcirc.S,1/2).
    % 
    % remarks:
    %   - The inverse is obtained by qcirc.T'
    T = qcirc.opFromMatrix([1 0; 0 (1+1i)/sqrt(2)],'T phase');

  end
  
  %% == Gate operator and projector getters ==
  
  methods(Static)

    % - Quantum states

    function v=b0
      % |0> qubit state.

      v = fockobj.bket('0');                                                    % Create ket
      v.name = '|0>';                                                           % Name it
    end

    function v=b1
      % |1> qubit state.

      v = fockobj.bket('1');                                                    % Create ket
      v.name = '|1>';                                                           % Name it
    end

    % - Physical gates of IBM Q
    
    function O=U1(lambda)
      % Gets the operator of IBM Q's physical U1 gate (1 qubit).
      %
      %   O = qcirc.U1(lambda)
      %
      % The operator matrix is
      %
      %   [1 0             ]
      %   [0 exp(1i*lambda)].
      %
      % arguments:
      %   lambda - phase angle in radians.
      %
      % returns:
      %   O      - The gate operator (a fockobj).
      %
      % See also fockobj, 
      % https://quantum-computing.ibm.com/support/guides/introduction-to-quantum-circuits

      om = [1 0; 0 exp(1i*lambda)];
      O  = qcirc.opFromMatrix(om,'U1 gate');
    end
    
    function O=U2(lambda,phi)
      % Gets the operator of IBM Q's physical U2 gate (1 qubit).
      %
      %   O = qcirc.U2(lambda,phi)
      %
      % The operator matrix is 
      %
      %   1/sqrt(2) * [1           -exp(1i*lambda)      ]
      %               [exp(1i*phi)  exp(1i*(lambda+phi))].
      %
      % arguments:
      %   lambda - phase angle #1 in radians.
      %   phi    - phase angle #2 in radians.
      %
      % returns:
      %   O      - The gate operator (a fockobj).
      %
      % See also fockobj,
      % https://quantum-computing.ibm.com/support/guides/introduction-to-quantum-circuits


      om = [1           -exp(1i*lambda)
            exp(1i*phi)  exp(1i*(lambda+phi))];
      O  = 1/sqrt(2)*qcirc.opFromMatrix(om,'U2 gate');
    end
    
    function O=U3(lambda,phi,theta)
      % Gets the operator of IBM Q's physical U3 gate (1 qubit).
      %
      %   O = qcirc.U3(lambda,phi,theta)
      %
      % The operator matrix is 
      %
      %   [ cos(theta/2)             -exp(1i*lambda)*sin(theta/2)       ]
      %   [ exp(1i*phi)*sin(theta/2)  exp(1i*lambda+1i*phi)*cos(theta/2)].
      %
      % arguments:
      %   lambda - phase angle #1 in radians.
      %   phi    - phase angle #2 in radians.
      %   theta  - phase angle #3 in radians.
      %
      % returns:
      %   O      - The gate operator (a fockobj).
      %
      % See also fockobj,
      % https://quantum-computing.ibm.com/support/guides/introduction-to-quantum-circuits

      om = [ cos(theta/2)             -exp(1i*lambda)*sin(theta/2)
             exp(1i*phi)*sin(theta/2)  exp(1i*lambda+1i*phi)*cos(theta/2)];
      O  = qcirc.opFromMatrix(om,'U3 gate');
    end
    
    function O=CNOT(varargin)
      % CNOT (controlled-NOT) gate operator (2 qubits).
      %
      %   O = qcirc.CNOT(s)
      %   O = qcirc.CNOT(N,c,t)
      %
      % arguments:
      %   s - Character array specifying the gate layout, a combination of
      %       - '*': Control qubit (exactly once),
      %       - 'U': Target qubit (exactly once), and
      %       - '-': Run-through qubit (arbitrary count).
      %       Qubits are ordered left to right, i.e., q[n]...q[0].
      %   N - Number of qubits in circuit.
      %   c - Zero-based index of control qubit.
      %   t - Zero-based index of target qubit.
      %
      % returns:
      %   O - The gate operator (a fockobj).
      %
      % See also fockobj

      % Intialization and checks                                                % -------------------------------------
      switch nargin                                                             % Branch for number of arguments >>
        case 1                                                                  %   One: -> first method signature
          s = varargin{1};                                                      %     Layout specifier
          if ~ischar(s) || isempty(s)                                           %     Layout spec. empty or not char.>>
            error('Bad layout, must be a non-empty character array.');          %       Error
          end                                                                   %     <<
          m = nnz(s=='*');                                                      %     Count controls
          if m==0; error('No control qubit specified in layout.'); end          %     No controls -> error
          if m>1; error('More than one control qubit specified in layout.'); end%     More than one control -> error
          m = nnz(s=='U');                                                      %     Count targets
          if m==0; error('No target qubit specified in layout.'); end           %     No targets -> error
          if m >1; error('More than one target qubit specified in layout.'); end%     More than one control -> error
          for i=1:length(s)                                                     %     Loop over chars. in layout spec. >>
            if ~ismember(s(i),'-*U')                                            %       Invalid character >>
              error('Invalid character ''%c'' in layout',s(i));                 %         Error
            end                                                                 %       <<
          end                                                                   %     <<
        case 3                                                                  %   Three: -> second method signature
          N = varargin{1};                                                      %     Number of qubits
          c = varargin{2};                                                      %     Control qubit
          t = varargin{3};                                                      %     Target qubit
          if (c<0 || c>=N); error('Bad control qubit index.'); end              %     Check control
          if (t<0 || t>=N); error('Bad target qubit index.'); end               %     Check target
          if (c==t)                                                             %     Control and target equal >>
            error('Control and target qubits cannot be identical.');            %       Error
          end                                                                   %     <<
          s      = char(ones(1,N)*'-');                                         %     Make layout specifier for qcirc.U
          s(N-c) = '*';                                                         %     ...
          s(N-t) = 'U';                                                         %     ...
        otherwise                                                               %   All other cases:
          error('Impermissible number of arguments.');                          %     Error
      end
      
      % Create operator                                                         % -------------------------------------
      O = qcirc.CU(qcirc.X,s);                                                  % Create controlled Pauli-X gate op.
      s = strrep(s,'*',char(8226));                                             % Have some fun with Unicode chars.
      s = strrep(s,'U',char(8853));                                             % ...
      O.name = sprintf('CNOT (%s)',s);                                          % Name gate operator

    end

    % - Fantasy gates I'd love to have physically

    function O=U(om,name)
      % Universal gate operator.
      % 
      %   O = qcirc.U(om)
      %   O = qcirc.U(om,name)
      %
      % arguments:
      %   om   - A unitary gate operator matrix in one of the forms
      %
      %          * 1 qubit gate:
      %
      %            om  |     <0|      <1|
      %            ----+-----------------
      %            |0> | om(1,1)  om(1,2)
      %            |1> | om(2,1)  om(2,2)
      %
      %          * 2 qubits gate:
      %
      %            om   |    <00|     <01|     <10|     <11|
      %            -----+-----------------------------------
      %            |00> | om(1,1)  om(1,2)  om(1,3)  om(1,4)
      %            |01> | om(2,1)  om(2,2)  om(2,3)  om(2,4)
      %            |10> | om(3,1)  om(3,2)  om(3,3)  om(3,4)
      %            |11> | om(4,1)  om(4,2)  om(4,3)  om(4,4)
      %
      %          * etc.
      %   name - Operator name (optional).
      %
      % returns:
      %   O    - The gate operator (a fockobj).
      %
      % throws exception:
      %   - If operator matrix is not square.
      %   - If operator matrix dimension is not a positive integer power of 2.
      %   - If operator matrix is not unitary.
      %
      % See also fockobj
      
      if nargin<2; name=''; end                                                 % Default name is empty
      try                                                                       % Try to >>
        O = qcirc.opFromMatrix(om,name);                                        %   Create operator from matrix
      catch e                                                                   % << Failed >>
        error(e.message);                                                       %   Error
      end                                                                       % <<
      try                                                                       % Try to >>
        assert(qcirc.isUnitary(O));                                             %   Check operator
      catch e                                                                   % << Failed >>
        error('Matrix is not unitary.');                                        %   Error
      end                                                                       % <<

    end
    
    function O=CU(U,s)
      % Controlled universal gate operator.
      %
      %   O = qcirc.CU(U,layout)
      %
      % arguments:
      %   U - Gate operator to control, single or multi qubit.
      %   s - Character array specifying the gate layout, a combination of
      %       - '*': Control qubit (exactly once),
      %       - 'U': Gate to control (exactly once), and
      %       - '-': Run-through qubit (arbitrary count).
      %       Qubits are ordered left to right, i.e., q[n]...q[0].
      %
      % returns:
      %   O - The gate operator (a fockobj).
      %
      % examples:
      %
      %   - CNOT gate with control q[0] and target q[1]:
      %
      %       CCNOT01 = qcirc.CU(qcirc.X,'U*')
      %
      %   - CNOT gate with control q[2], target q[0], and run-through q[1]:
      %
      %       CCNOT20 = qcirc.CU(qcirc.X,'*-U')
      %
      %   - Calls may be nested. E.g., a Toffoli gate with controls q[1], q[2] 
      %     and target q[0] is obtained by:
      %
      %       CCCNOT120 = qcirc.CU(qcirc.CU(qcirc.X,'*U'),'*U')
      %
      % See also fockobj

      % Initialization and checks                                               % -------------------------------------
      if ~ischar(s) || isempty(s)                                               % Layout specifier empty or not char.>>
        error('Bad layout, must be a non-empty character array.');              %   Error
      end                                                                       % <<
      n = length(s);                                                            % Get length of layout specifier
      m = nnz(s=='*');                                                          % Count controls
      if m==0; error('No control qubits specified in layout.'); end             % No controls -> error
      if m >1; error('More than one control qubit specified in layout.'); end   % More than one control -> error
      if ~ismember('U',s); warning('No target specified in layout.'); end       % No target(s) -> warning
      for i=1:n                                                                 % Loop over chars. in layout spec. >>
        if ~ismember(s(i),'-*U')                                                %   Invalid character >>
          error('Invalid character ''%c'' in layout',s(i));                     %     Error
        end                                                                     %   <<
      end                                                                       % <<
      [~,nU] = fockbasis(U).getNumSectors();                                    % Get number of qubits in gate U
      idU = qcirc.id; for i=2:nU; idU = kron(idU,qcirc.id); end                 % Make nU-qubit identity operator
      O1  = 1;                                                                  % Left summand of operator
      O2  = 1;                                                                  % Right summand of operator

      % Create operator                                                         % -------------------------------------
      for i=0:n-1                                                               % Loop over circuit qubits >>
        if s(n-i)=='*'                                                          %   q[n-i] is the control qubit >> 1)
          O1 = kron(qcirc.b0*qcirc.b0',O1);                                     %     Build left summand of operator
          O2 = kron(qcirc.b1*qcirc.b1',O2);                                     %     Build right summand of operator
        elseif s(n-i)=='U'                                                      %   << q[n-1] is the target qubit >> 1)
          O1 = kron(idU,O1);                                                    %     Build left summand of operator
          O2 = kron(  U,O2);                                                    %     Build right summand of operator
        else                                                                    %   << Other qubits >>
          O1 = kron(qcirc.id,O1);                                               %     Build left summand of operator
          O2 = kron(qcirc.id,O2);                                               %     Build right summand of operator
        end                                                                     %   <<
      end                                                                       % <<
      O = O1+O2;                                                                % Create operator
      O.name = strrep(string(s),'U',sprintf('(%s)',U.name));                    % Name operator
      
      % Remarks
      % 1) (N-i) mirrors the order of the tensor product factors as q[N] is
      %    the _first_ and q[0] is the _last_ factor.
    end

    function O=BOOL(tt,layout)
      % Universal Boolean gate operator.
      %
      %   O = qcirc.BOOL(tt,layout)
      %
      % arguments:
      %   tt - Truth table, a logical vector of 2^N elements. N denotes the 
      %        number of arguments of the Boolean function.
      %   s  - Character vector specifying the gate layout, a combination of
      %        - 'i': Input qubit (exactly N times),
      %        - 'o': Output qubit (exactly once), and
      %        - '-': Run-through qubit (arbitrary count).
      %        Qubits are ordered left to right, i.e., q[n]...q[0].
      %
      % returns:
      %   O - The gate operator (a fockobj).
      %
      % examples:
      %
      %   - Toffoli gate operator: q[2] = q[1] AND q[0]
      % 
      %       tt = [ false false false true ];
      %       O = qcirc.booleanOp(tt,'oii');
      %
      % See also fockobj
      
      % Check and prepare                                                       % -------------------------------------
      % - Layout specifier                                                      % 
      if ~isvector(layout) || ~ischar(layout)                                   % Check type of layout arg. >>
        error('Gate layout specifier must be a character vector.');             %   Invalid -> error
      end                                                                       % <<
      iid = find(layout=='i');                                                  % Input qubit indices in layout spec.
      oid = find(layout=='o');                                                  % Output qubit indices in layout spec.
      rid = find(layout=='-');                                                  % Run-thru qubit indices in layout sp.
      if length(oid) ~= 1                                                       % More or less than one output qubit >> 
        error('Gate layout specifier must contain exactly one ''o''.');         %   Error
      end                                                                       % <<
      if length(iid)+length(oid)+length(rid) ~= length(layout)                  % Invalid characters in layout >>
        error('Gate layout specifier must ''o'', ''i'', and ''-'' only.');      %   Error
      end                                                                       % <<
      
      % - Thruth table                                                          % -------------------------------------
      if ~isvector(tt) || ~islogical(tt)                                        % Check type of tt argument >>
        error('Thruth table must be a logical vector.');                        %   Invalid -> error
      end                                                                       % <<
      Ntt = log2(length(tt));                                                   % # input qubits to function from tt
      if Ntt ~= round(Ntt)                                                      % Not an integer power of 2 >>
        error('Length of truth table vector must be an integer power of 2.');   %   Error
      end                                                                       % <<
      
      % - Get and check numbers of qubits                                       % -------------------------------------
      Nq = length(layout);                                                      % # qubits in gate
      Ni = length(iid);                                                         % # input qubits to funct. from layout
      if Ni ~= Ntt                                                              % Inconsistent number of inputs >>
        error('Length of truth table vector does not match layout specifier.'); %   Error
      end                                                                       % <<
      
      % Make map of input bit combinations to output values of Boolean function % -------------------------------------    
      ttmap = containers.Map();                                                 % Map of bit combinations to truth val.
      for i=1:length(tt)                                                        % Loop over truth table entries >>
        key = qcirc.int2binstr(i-1,Ni);                                         %   Get bit combination of inputs
        val = tt(i);                                                            %   Get value from truth table
        ttmap(key) = val;                                                       %   Store key-value pair
      end                                                                       % <<
      
      % Build operator element-wise                                             % ------------------------------------- 
      O    = 0;                                                                 % Zero-initialize operator
      tdel = repmat('�',1,Nq);                                                  % Aux. array of tensor prod. delimiters
      for i=0:2^Nq-1                                                            % Loop over bit combinations >>
        bc  = qcirc.int2binstr(i,Nq);                                           %   Get bit combination as char. array
        in  = bc(iid);                                                          %   Get input bit combination
        out = bc(oid);                                                          %   Get output bit
        f   = ttmap(in);                                                        %   Get Boolean fnc. value for input
        if f                                                                    %   Boolean function value is true >>
          if out=='0'; out='1'; else; out='0'; end                              %     Flip output bit
        end                                                                     %   <<
        ketid = bc;                                                             %   Ket <- bit combination
        braid = bc; braid(oid)=out;                                             %   Bra <- bit comb. w/ modified output
        ketid = reshape([ketid; tdel],1,2*Nq); ketid = ketid(1:2*Nq-1);         %   Make basis vector identifier of ket
        braid = reshape([braid; tdel],1,2*Nq); braid = braid(1:2*Nq-1);         %   Make basis vector identifier of bra
        ket = fockobj.bket(ketid);                                              %   Create ket
        bra = fockobj.bket(braid)';                                             %   Create bra
        O = O + ket*bra;                                                        %   Contribute element to operator
      end                                                                       % <<
      
      % Final checks                                                            % -------------------------------------
      assert(qcirc.isUnitary(O)); % Remove later                                % Check newly created operator
      
    end
    
    % - Sub-space projectors
    
    function P=proj(s)
      % Creates projector on a sub-space of a quantum cicuit's state space.
      %
      %   P = qcirc.proj(spec)
      %
      % arguments:
      %   s - Character array specifying the sub-space, a combination of
      %       - '0': Qubit is in state |0>,
      %       - '1': Qubit is in state |1>, and
      %       - '*': Qubit is in any state.
      %       Qubits are ordered left to right, i.e., q[n]...q[0].
      %
      % returns:
      %   P - The projector (a fockobj).
      % 
      % example:
      %
      %   qcirc.proj('**0*1') returns the projector of the sub-space of the
      %   state space of a five qubit circuit where q[0]=|1> and q[2]=|0>.

      N = strlength(s);                                                         % Get number of qubits
      P = 1;                                                                    % Initialize projector
      for i=0:N-1                                                               % Loop over qubits >>
        if     s(N-i)=='0'; P = kron(qcirc.b0*qcirc.b0',P);                     %   Project on q[i]=|0>
        elseif s(N-i)=='1'; P = kron(qcirc.b1*qcirc.b1',P);                     %   Project on q[i]=|1>
        else              ; P = kron(qcirc.id          ,P);                     %   Do not project on q[i]
        end                                                                     %   ...
      end                                                                       % <<
      P.name = sprintf('Sub-space projector (%s)',s);                           % Name projector

    end

  end

  %% == Measurement ==
  
  methods(Static)

    function probe(psi,varargin)
      % Probes the state of qubits in a quantum circuit and displays the result 
      % in a histogram.
      %
      %   qcirc.probe(psi)
      %   qcirc.probe(___,Name,Value)
      %
      % arguments:
      %   Psi      - State of quantum circuit, a fockobk (ket).
      %   ___      - Any other argument list.
      % 
      % returns:
      %   nothing
      %
      % additional properties:
      %   'qubits' - A vector of zero-based indexes of the qubits to probe,
      %              default is all qubits.
      %
      %   'shots'  - Number of simulated shots, default is 1024.
      %
      %   'name'   - A name (displayed as histogram title).
      
      % Initialize and check                                                    % -------------------------------------
      if ~isa(psi,'fockobj') || psi.getType()~=fock.OBJ_KET                     % Psi not a ket vector >>
        error('Psi must be ket fockobj.');                                      %   Error
      end                                                                       % <<
      [~,N] = fockbasis(psi).getNumSectors();                                   % Get number of qubits in circuit
      try qcirc.checkOptions(varargin,{'qubits','shots','name'});               % Check option names
      catch e; error(e.message); end                                            % Invalid option name(s) -> error
      qubits = qcirc.getOption(varargin,'qubits',N-1:-1:0);                     % Get indexes of qubits to probe
      if ~all(ismember(qubits,0:N-1))                                           % There are invalid qubit indexes >>
        error('Qubit indexes must be integers from interval 0:%d.',N-1);        %   Error
      end                                                                       % <<
      shots = qcirc.getOption(varargin,'shots',1024);                           % Get number of simulated shots
      if ~isnumeric(shots) || ~isscalar(shots) || rem(shots,1)~=0 || shots<1    % Invalid number of shots >>
        error('''shots'' must be a positive integer.');                         %   Error
      end                                                                       % <<
      NQ = length(qubits);                                                      % Get number of probe qubits
      name = qcirc.getOption(varargin,'name',psi.name);                         % Get name
      
      % Compute theoretical result                                              % -------------------------------------
      NN = 2^NQ;                                                                % Number of qubit sub-spaces
      cats = cell(1,NN);                                                        % Histogram category names
      probs = zeros(1,NN);                                                      % Theoretical probabilities
      for i=0:NN-1                                                              % Loop over probe qubit sub-spaces >>
        cbits = qcirc.int2binstr(i,NQ);                                         %   Get probe qubit combination
        cats{i+1} = cbits;                                                      %   Set histogram category name
        pspec = char(ones(1,N)*'*');                                            %   Make sub-space projector specifier
        for j=0:NQ-1; pspec(N-qubits(NQ-j))=cbits(NQ-j); end                    %   ...
        P = qcirc.proj(pspec);                                                  %   Make sub-space projector
        probs(i+1) = psi'*P*psi;                                                %   Compute sub-space probability
      end                                                                       % <<

      % Simulate shot statistics                                                % -------------------------------------
      counts = zeros(1,NN);                                                     % Array of bin counts
      limits = zeros(1,NN+1);                                                   % Array of probability bin limits ...
      for i=1:NN; limits(1,i+1) = limits(1,i+0)+probs(1,i); end                 % ... Initialize
      for i=1:shots                                                             % Simulate shots >>
        x = rand();                                                             %   Get a random "probability" value
        j = find(circshift(x<limits,-1).*(x>=limits));                          %   Find probability bin index
        counts(1,j) = counts(1,j)+1;                                            %   Increment count in this bin
      end                                                                       % <<

      % Show histgrams                                                          % -------------------------------------
      probs   = round(probs*NN*shots);                                          % Theoretical probabilities -> "counts"
      options = { 'Categories'   , cats         , ...                           % Histogram options
                  'Normalization', 'probability', ...
                };                                                              % ...
      histogram('BinCounts',counts,options{1:end},'BarWidth',0.8);              % Show histogram of simulated counts
      hold on                                                                   % Stay in plot
      histogram('BinCounts',probs,options{1:end},'BarWidth',0.4);               % Add theoretical probabilities
      legend(sprintf('simulation of %d shots',shots), ...                       % Set histogram legend
        'theoretical probabilities');                                           % ...
      s = 'Quantum circuit state';                                              % Make figure title
      if ~isempty(name); s = s + " """ + string(name) + """"; end               % ...
      title(s);                                                                 % Set axes title
      hold off                                                                  % Stop staying in plot
      
    end

  end
  
  %% == Auxiliary methods ==
  
  methods(Static)

    function TF=isUnitary(O)
      % Checks if an operator is unitary.
      %
      %   TF=qcirc.isUnitary(O)
      %
      % arguments:
      %   O  - The operator, a fockobj.
      %
      % returns:
      %   TF - true if the operator matrix is unitary, false otherwise.
      %
      % See also fockobj
      
      om = fockbasis(O).realize(O);                                             % Get operator matrix
      om(1,:) = []; om(:,1) = [];                                               % Remove vacuum row and column
      TF = all(all(om'*om-eye(length(om))<eps));                                % Test if unitary
    end

    function disp(O)
      % Pretty-prints a quantum gate operator or a quantum state.
      %
      %   qcirc.disp(O)
      %   qcirc.disp(O,prec)
      %
      % arguments:
      %   O    - The operator or quantum state, a fockobj.
      %   prec - Number of significant digits to print (optional, default is 4).
      %
      % returns:
      %   nothing
      %
      % See also fockobj

      if nargin<2; prec=4; end                                                  % Default number of sign. digits
      [~,N] = fockbasis(O).getNumSectors();                                     % Get number of qubits
      switch O.getType()                                                        % Branch for fockobj type >>
        case fock.OBJ_LOP; s = 'quantum gate operator';                         %   Operator
        case fock.OBJ_KET; s = 'quantum state';                                 %   Ket
        case fock.OBJ_BRA; s = 'quantum state (bra)';                           %   Bra
        otherwise; error('Invalid argument.');                                  %   Otherwise -> error
      end                                                                       % <<
      fprintf('  %s <strong>%s</strong>\n\n',fock.HREF('qcirc'),s);             % Print headline
      if ~strcmp(string(O.name),"")                                             % Non-empty name >>
        fprintf('    Name  : <strong>%s</strong>\n',O.name);                    %   Print name
      end                                                                       % <<
      fprintf('    Qubits: <strong>%d</strong>\n\n',N);                         % Print number of qubits
      fprintf(qcirc.sprintqo(O,'prec',prec));                                   % Print operator matrix
      fprintf('\n');                                                            % Print line break
    end

    function latex(O,prec)
      % Generates LaTeX code for pretty-printing a quantum gate operator or a 
      % quantum state.
      %
      %   qcirc.latex(O)
      %   qcirc.latex(O,prec)
      %
      % arguments:
      %   O    - The operator or quantum state, a fockobj.
      %   prec - Number of significant digits to print (optional, default is 4).
      %
      % returns:
      %   nothing
      %
      % See also fockobj
      
      if nargin<2; prec=4; end                                                  % Default number of sign. digits
      fprintf('%% Generated by FockBox'' qcirc.latex method\n');                % LaTeX preamble
      fprintf(qcirc.sprintqo(O,'prec',prec,'mode','latex'));                    % Print operator matrix
      fprintf('\n');                                                            % Print line break
    end

    function s=sprintqo(O,varargin)
      % Pretty-prints a quantum gate operator or a quantum state into a string.
      % 
      %   s = qcirc.sprintqo(O)
      %   s = qcirc.sprintqo(___,Name,Value)
      %
      % arguments:
      %   O           - The operator or quantum state, a fockobj.
      %   prec        - Number of significant digits to print (optional, default 
      %                 is 4).
      %   Name,Value  - Additional properties for printing.
      %   ___         - Any other argument list.
      %
      % returns:
      %   s           - The string representation of the quantum gate
      %                 operator.
      %
      % additional properties:
      %   'prec'      - Number of significant digits to print. The default
      %                 value is 4.
      % 
      %   'mode'      - Printing mode, value can be one of the following
      %     'console' - For Matlab console (default)
      %     'latex'   - LaTeX code
      %     'plain'   - Plain ASCII
      %
      % See also fockobj

      % Initialize                                                              % -------------------------------------
      try qcirc.checkOptions(varargin,{'prec','mode'});                         % Check option names
      catch e; error(e.message); end                                            % Invalid option name(s) -> error
      mode = qcirc.getOption(varargin,'mode','console');                        % Get printing mode
      if ~ismember(mode,{'console'; 'latex'; 'plain'})                          % Invalid mode value >>
        error('''mode'' must be one of ''console'', ''latex'', or ''plain''.'); %   Error
      end                                                                       % <<
      prec = qcirc.getOption(varargin,'prec',4);                                % Get no. of significant digits
      if ~isnumeric(prec) || ~isscalar(prec) || rem(prec,1)~=0 || prec<0        % Invalid precision value >>
        error('''prec'' must be a non-negative integer');                       %   Error
      end                                                                       % <<
      s     = "";                                                               % Initialize output string
      B     = fockbasis(O);                                                     % Get operator basis
      [~,N] = B.getNumSectors();                                                % Get number of qubits
      om    = qcirc.opToMatrix(O);                                              % Get operator matrix

      % Preflight: Determine column widths                                      % -------------------------------------
      cw = zeros(size(om,2)+1,1);                                               % Initialize column with array
      for j=1:size(om,1)                                                        % Loop over op. matrix columns >>
        [~,l]   = qcirc.sprintvec(qcirc.int2binstr(j-1,N),'bra',mode);          %   Get printable length
        cw(1)   = max(cw(1),l);                                                 %   Aggregate with of header column
        cw(j+1) = l;                                                            %   Initialize with if header row
      end                                                                       % <<
      for i=1:size(om,1)                                                        % Loop over op. matrix rows >>
        for j=1:size(om,2)                                                      %   Loop over op. matrix columns >>
          [~,l] = qcirc.sprintcmplx(om(i,j),prec,mode);                         %     Get printable len. of val. string
          cw(j+1) = max(cw(i+1),l);                                             %     Aggregate column width
        end                                                                     %   <<
      end                                                                       % <<

      % Pretty-print operator matrix                                            % -------------------------------------
      if strcmp(mode,'latex')                                                   % Latex mode >>
        sRH = '    '; sCD = ' & '; sRT = ' \\\\\n';                             %   Set table delimiters
        s = s + '%% \\usepackage{rotating}\n';                                  %   ...
        s = s + '{\n';                                                          %   ...
        s = s + '  \\def\\s{\\scriptstyle}\n';                                  %   ...
        s = s + '  \\def\\b#1{\\mathbf{#1}}\n';                                 %   ...
        s = s + '  \\def\\e{\\mathrm{e}}\n';                                    %   ...
        s = s + '  \\def\\i{\\mathrm{i}}\n';                                    %   ...
        s = s + '  \\def\\ket#1{\\s|#1\\rangle}\n';                             %   ...
        s = s + '  \\def\\bra#1{\\rotatebox{90}{$\\s\\langle #1|$}}\n';         %   ...
        s = s + '  \\hspace*{-0.4ex}\\left(\\hspace*{-0.4ex}\\begin{array}{';   %   ...
        if size(om,1)>1; s = s + 'c|'; end                                      %   ...
        for j=1:size(om,2); s = s + 'c'; end                                    %   ...
        s = s + '}\n';                                                          %   ...
      else                                                                      % << All other modes >>
        sRH = '    '; sCD = ' '; sRT = '\n';                                    %   Set table delimiters
      end                                                                       % <<
      if size(om,2)>1                                                           % Not a ket >>
        s = s + sRH;                                                            %   Indent line
        if size(om,1)>1                                                         %   Not a bra >>
          f = '%'+sprintf("%ds",cw(1)); s = s  + sprintf(f,'');                 %     Output header of first row
        end                                                                     %   <<
        for j=1:size(om,2)                                                      %   Loop over op. matrix columns >>
          t = qcirc.sprintvec(qcirc.int2binstr(j-1,N),'bra',mode,cw(j+1));      %     Print bra vector
          if j>1 || size(om,1)>1; s = s + sCD; end                              %     Output column delimiter
          s = s + t;                                                            %     Output column header
        end                                                                     %   <<
        if strcmp(mode,'latex'); s = s + ' \\\\\\hline\n'; else; s = s +sRT; end%   Output line break
      end                                                                       % <<
      for i=1:size(om,1)                                                        % Loop over op. matrix rows >>
        s = s + sRH;                                                            %   Indent line
        if size(om,1)>1                                                         %   Not a bra >>
          s = s + qcirc.sprintvec(qcirc.int2binstr(i-1,N),'ket',mode,cw(1));    %     Output row header
        end                                                                     %   <<
        for j=1:size(om,2)                                                      %   Loop over op. matrix columns >>
          if j>1 || size(om,1)>1; s = s + sCD; end                              %     Output column delimiter
          s = s + qcirc.sprintcmplx(om(i,j),prec,mode,cw(j+1));                 %     Output values
        end                                                                     %   <<
        if i<size(om,1); s = s + sRT; end                                       %   Output line break (except last row)
      end                                                                       % <<
      if strcmp(mode,'latex')                                                   % Latex mode >>
        s = s + '\n  \\end{array}\\!\\right)\n}';                               %   LaTeX postamble
      end                                                                       % <<
    end
    
  end
  
  %% == Auxiliary private methods ==

  methods(Static,Access=private)
    
    function O=opFromMatrix(om,name)
      % Builds an operator from an operator matrix.
      % 
      %   O = qcirc.opFromMatrix(om)
      %   O = qcirc.opFromMatrix(om,name)
      %
      % arguments:
      %   om   - A 2^n by 2^n operator matrix with n>0.
      %   name - Operator name (optional).
      %
      % returns:
      %   O    - The gate operator (a fockobj).
      %
      % throws exception:
      %   - If operator matrix is not square.
      %   - If operator matrix dimension is not a positive integer power of 2.
      %
      % See also opToMatrix, fockobj

      % Checks and initialization                                               % -------------------------------------
      assert(size(om,1)==size(om,2),'Matrix mut be square.');                   % Check matrix is square
      N = log2(size(om,1));                                                     % Get number of qubits
      assert(mod(N,1)==0&&N>0,...                                               % Check matrix dimension is 2^N, N>0
        'Matrix dimension must be a positive integer power of 2.');             % ...
      
      % Make operator                                                           % -------------------------------------
      B  = sector(qcirc.basis^N,N);                                             % Get Fock space basis of N qubits
      om = [ zeros(1,2^N); om];                                                 % Insert heading row of zeros into om
      om = [ zeros(2^N+1,1) om];                                                % Insert heading col.of zeros into om
      O  = B.unrealize(om);                                                     % Convert matrix to Fock space operator
      if nargin>1; O.name = name; end                                           % Set operator name

      % Old code for 2x2 matrices                                               % -------------------------------------
      %b0 = qcirc.b0;                                                           % Get |0>
      %b1 = qcirc.b1;                                                           % Get |1>
      %O  = om(1,1)*(b0*b0') + om(1,2)*(b0*b1')...                              % Assemble operator point-wise
      %   + om(2,1)*(b1*b0') + om(2,2)*(b1*b1');                                % ...
    end

    function om=opToMatrix(O)
      % Gets the operator matrix of a gate operator.
      % 
      %   om = qcirc.opToMatrix(O)
      %
      % arguments:
      %   O  - The gate operator,a fockobj.
      %
      % returns:
      %   om - The operator matrix. See documentation of opFromMatrix for
      %        format.
      %
      % See also opFromMatrix, fockobj

      [~,N] = fockbasis(O).getNumSectors();                                     % Get number of qubits
      B  = sector(qcirc.basis^N,N);                                             % Get Basis of Fock space sector
      om = full(B.realize(O));                                                  % Get operator matrix
      switch O.getType()                                                        % Branch for fockobj type >>
        case fock.OBJ_LOP; om = om(2:end,2:end);                                %   Operator
        case fock.OBJ_BRA; om = om(1:end,2:end);                                %   Bra (undocumented feature)
        case fock.OBJ_KET; om = om(2:end,1:end);                                %   Ket (undocumented feature)
        otherwise; error('Argument is no Fock space operator.');                %   Otherwise -> error
      end                                                                       % <<

    end

    function [s,len]=sprintcmplx(v,prec,mode,width)
      % Thrifty-prints a complex value to a formatted string (with markup).
      %
      %   [s,len] = sprintcmplx(v,prec,mode,width)
      %
      % arguments:
      %   v     - The value to print.
      %   prec  - Number of significant digits to print.
      %   mode  - Mode, one of 'console', 'latex', or 'plain'.
      %   width - With of printed string (padded with spaces).
      %
      % returns:
      %   s     - The formatted string
      %   len   - Number of printable characters in string.

      % Print complex number                                                    % -------------------------------------
      f  = "%"+sprintf('.%dg',prec);                                            % Format string w/ precision
      fs = "%"+sprintf('+.%dg',prec);                                           % Format string w/ precision and sign
      if isreal(v)                                                              % Value is real >>
        s = sprintf(f,real(v));                                                 %   Print real part only
      elseif isreal(v*1i)                                                       % << Value is imaginary >>
        s = sprintf(f+'i',imag(v));                                             %   Print imaginary part only
      else                                                                      % << Value is complex >>
        s = sprintf(f+fs+'i',real(v),imag(v));                                  %   Print complex value
      end                                                                       % <<
      len = strlength(s);                                                       % Get # printable characters
      if nargin<4; width=len; end                                               % Minimal print with

      % Add mode specific markup                                                % -------------------------------------
      f = "%"+sprintf('%ds',width);                                             % Make print format
      if strcmp(mode,'latex')                                                   % LaTeX mode >>
        s = strrep(s,'i','\i');                                                 %   Replace imaginary unit
        if ~strcmp(s,"0"); s = '\b{'+s+'}'; end                                 %   Print non-zero value bold
        s = sprintf(f,s);                                                       %   Space-pad printed value
        len = strlength(s);                                                     %   Everything is protable
        s = replace(s,'\','\\');                                                %   Escape backslash
      else                                                                      % << Matlab console or plain mode >>
        if ~strcmp(s,"0") && strcmp(mode,'console')                             %   Value not zero and console mode >>
          s = '<strong>'+sprintf(f,s)+'</strong>';                              %     Print value in bold face
        else                                                                    %   << Value zero >>
          s = sprintf(f,s);                                                     %     Print value in normal face
        end                                                                     %   <<
      end                                                                       % <<

    end

    function [s,len]=sprintvec(name,type,mode,width)
      % Prints a ket or bra identifier to a formatted string (with markup).
      %
      %   [s,len] = sprintvec(name,type,mode,width)
      %
      % arguments:
      %   name  - Name of the vector.
      %   type  - Vector type, 'ket' or 'bra'.
      %   mode  - Mode, one of 'console', 'latex', or 'plain'.
      %   width - With of printed string (padded with spaces).
      %
      % returns:
      %   s     - The formatted string
      %   len   - Number of printable characters in string.
      
      if strcmp(type,'bra')                                                     % Print a bra vector >>
        sl = '<'; sr = '|';                                                     %   Set delimiters
        if strcmp(mode,'latex'); sl = '\bra{'; sr='}'; end                      %   Set LaTeX left delimiter
      else                                                                      % << Print a ket vector >>
        sl = '|'; sr = '>';                                                     %   Set delimiters
        if strcmp(mode,'latex'); sl = '\ket{'; sr='}'; end                      %   Set LaTeX right delimiter
      end                                                                       % <<
      if nargin<4; width=strlength(name); end                                   % Minimal print with
      f = "%"+sprintf('%ds',width);                                             % Make print format
      s = sprintf(f,string(sl)+string(name)+string(sr));                        % Print ket
      len = strlength(s);                                                       % Printable length
      s = replace(s,'\','\\');                                                  % Escape backslash

    end

    function s=int2binstr(n,len)
      % Converts an non-negative integer to a binary number string.
      %
      %   s = qcirc.int2binstr(n)
      %   s = qcirc.int2binstr(n,len)
      %
      % arguments:
      %   n   - The integer.
      %   len - The minimal length of the binary number (padded with
      %         leading '0's)
      %
      % returns:
      %   s   - The binary number represented by a character array.

      if nargin<2; len=1; end                                                   % Default minimal length
      s = "";                                                                   % Initialize result
      while n>0                                                                 % Value (still) greater than zero >>
        b = rem(n,2);                                                           %   Remainder of dividing by 2
        s = string(b) + s;                                                      %   Prefix result with rem. as char.
        n = floor(n/2);                                                         %   Integer division by 2
      end                                                                       % <<
      while true                                                                % Zero padding loop >>
        if strlength(s)>=len; break; end                                        %   Result long enough -> break
        s = "0"+s;                                                              %   Prefix result by "0"
      end                                                                       % <<
      s = convertStringsToChars(s);                                             % Convert result to character array

    end

    function v=getOption(nvlist,name,default)
      % Gets an option value from a name-value list.
      %
      %   v = qcirc.getOption(name,default,list)
      %
      % arguments:
      %   nvlist  - The name-value list, a cell vector.
      %   name    - The option name.
      %   default - The default value.
      %
      % returns:
      %   v       - The value, default if name is not in the list.

      v = default;                                                              % Set result to default
      for i=1:2:length(nvlist)                                                  % Loop over name-value pairs >>
        if strcmp(nvlist{i},name)                                               %   Name found >>
          v = nvlist{i+1};                                                      %     Get value
          return;                                                               %     That was it...
        end                                                                     %   <<
      end                                                                       % <<

    end
    
    function checkOptions(nvlist,nlist)
      % Checks names in an option list.
      %
      %   qcirc.checkOptions(nvlist,nlist));
      %
      % arguments:
      %   nvlist - A name-value list.
      %   nlist  - A list if permissible names.
      %
      % returns:
      %   nothing
      %
      % throws exception:
      %   - If name-value list containes names not mentioned in name list.
      
      for i=1:2:length(nvlist)                                                  % Loop over name-value pairs >>
        assert(ismember(nvlist{i},nlist), ...                                   %   Assert name is in list
          'Unrecognized option ''%s''.',nvlist{i});                             %   ...
      end                                                                       % <<
    end
    
  end

  %% == Strange, draft or dubious... ==

  methods(Static)
    
    % Nothing in this category presently
    
  end

end

% EOF