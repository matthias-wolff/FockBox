classdef ibmqx
  % Library of IBM Quantum Experience's quantum gates as Fock space
  % operators.
  %
  % TODO:
  % - Implement measurement!
  
  %% == Constants ==
  
  properties (Constant)

    % |0> qubit state.
    b0 = fockobj.bket('0');

    % |1> qubit state.
    b1 = fockobj.bket('1');

    % The qubit basis.
    basis = fockbasis({ibmqx.b0 ibmqx.b1});
    
    % Identity gate operator (1 qubit).
    id = ibmqx.opFromMatrix([1 0; 0 1],'Identity gate');

    % Pauli-X gate operator (1 qubit).
    X = ibmqx.opFromMatrix([0 1; 1 0],'Pauli-X gate');

    % Pauli-Y gate operator (1 qubit).
    Y = ibmqx.opFromMatrix([0 -1i; 1i 0],'Pauli-Y gate');
  
    % Pauli-Z gate operator (1 qubit).
    Z = ibmqx.opFromMatrix([1 0; 0 -1],'Pauli-Z gate');

    % Hadamard gate operator (1 qubit).
    H = 1/sqrt(2)*ibmqx.opFromMatrix([1 1; 1 -1],'Hadamard gate');

    % S phase gate operator (1 qubit), S = mpowerf(ibmqx.Z,1/2).
    S = ibmqx.opFromMatrix([1 0; 0 1i],'S phase gate');

    % T phase gate operator (1 qubit), T = mpowerf(ibmqx.S,1/2).
    T = ibmqx.opFromMatrix([1 0; 0 (1+1i)/sqrt(2)],'T phase gate');

  end
  
  %% == Static API ==
  
  methods(Static)

    function O=opFromMatrix(om,name)
      % Builds a gate operator from an operator matrix.
      % 
      %   O = opFromMatrix(om)
      %   O = opFromMatrix(om,name)
      %
      % arguments:
      %   om   - The operator matrix in one of the forms
      %
      %          * 1 qubit operator:
      %
      %            om  | |0>      |1>
      %            ----+-----------------
      %            |0> | om(1,1)  om(1,2)
      %            |1> | om(2,1)  om(2,2)
      %
      %          * 2 qubits operator:
      %
      %            om   | |00>     |01>     |10>     |11>
      %            -----+-----------------------------------
      %            |00> | om(1,1)  om(1,2)  om(1,3)  om(1,4)
      %            |01> | om(2,1)  om(2,2)  om(2,3)  om(2,4)
      %            |10> | om(3,1)  om(3,2)  om(3,3)  om(3,4)
      %            |11> | om(4,1)  om(4,2)  om(4,3)  om(4,4)
      %
      %          * etc.
      %   name - Operator name (optionl)
      %
      % returns:
      %   O    - The gate operator (a fockobj).
      %
      % See also opToMatrix

      % Checks and initialization                                               % -------------------------------------
      assert(size(om,1)==size(om,2),'Matrix mut be square.');                   % Check matrix is square
      N = log2(size(om,1));                                                     % Get number of qubits
      assert(mod(N,1)==0&&N>0,...                                               % Check matrix dimension is 2^N, N>0
        'Matrix dimension must be a positive integer power of 2.');             % ...
      
      % Make operator                                                           % -------------------------------------
      B  = sector(ibmqx.basis^N,N);                                             % Get Fock space basis of N qubits
      om = [ zeros(1,2^N); om];                                                 % Insert heading row of zeros into om
      om = [ zeros(2^N+1,1) om];                                                % Insert heading col.of zeros into om
      O  = B.unrealize(om);                                                     % Convert matrix to Fock space operator
      if nargin>1; O.name = name; end                                           % Set operator name

      % Old code for 2x2 matrices                                               % -------------------------------------
      %b0 = ibmqx.b0;                                                           % Get |0>
      %b1 = ibmqx.b1;                                                           % Get |1>
      %O  = om(1,1)*(b0*b0') + om(1,2)*(b0*b1')...                              % Assemble operator point-wise
      %   + om(2,1)*(b1*b0') + om(2,2)*(b1*b1');                                % ...
    end

    function om=opToMatrix(O)
      % Gets the operator matrix of a gate operator.
      % 
      %   om = opToMatrix(O)
      %
      % arguments:
      %   O  - The gate operator (a fockobj).
      %
      % returns:
      %   om - The operator matrix. See documentation of opFromMatrix for
      %        format.
      %
      % See also opFromMatrix

      B  = fockbasis(O);
      om = full(B.realize(O));
      om = om(2:end,2:end);
    end
    
    function O=U1(lambda)
      % Gets the operator of IBM QX's first physical gate (1 qubit).
      %
      %   O = U1(lambda)
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

      om = [1 0; 0 exp(1i*lambda)];
      O  = ibmqx.opFromMatrix(om,'U1 gate');
    end
    
    function O=U2(lambda,phi)
      % Gets the operator of IBM QX's second physical gate (1 qubit).
      %
      %   O = U2(lambda,phi)
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

      om = [1           -exp(1i*lambda)
            exp(1i*phi)  exp(1i*(lambda+phi))];
      O  = 1/sqrt(2)*ibmqx.opFromMatrix(om,'U2 gate');
    end
    
    function O=U3(lambda,phi,theta)
      % Gets the operator of IBM QX's third physical gate (1 qubit).
      %
      %   O = U3(lambda,phi,theta)
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

      om = [ cos(theta/2)             -exp(1i*lambda)*sin(theta/2)
             exp(1i*phi)*sin(theta/2)  exp(1i*lambda+1i*phi)*cos(theta/2)];
      O  = ibmqx.opFromMatrix(om,'U3 gate');
    end
    
    function O=CNOT(N,c,t)
      % CNOT (controlled-NOT) gate operator (2 qubits).
      %
      %   O = CNOT(c,t)
      %
      % arguments:
      %   N - Number of qubits in circuit.
      %   c - Zero-based index of control qubit.
      %   t - Zero-based index of target qubit.
      %
      % returns:
      %   O - The gate operator (a fockobj).
      %
      % TODO:
      % - Implement more elegantly!
      
      % Checks
      if (c<0 || c>=N); error('Bad control qubit index.'); end
      if (t<0 || t>=N); error('Bad target qubit index.'); end
      if (c==t); error('Control and target qubits cannot be identical.'); end

      % Create operator
      O1 = 1;
      O2 = 1;
      for i=N-1:-1:0
        if i==c
          O1 = kron(O1,ibmqx.b0*ibmqx.b0');
          O2 = kron(O2,ibmqx.b1*ibmqx.b1');
        elseif i==t
          O1 = kron(O1,ibmqx.id);
          O2 = kron(O2,ibmqx.X );
        else
          O1 = kron(O1,ibmqx.id);
          O2 = kron(O2,ibmqx.id);
        end  
      end
      O = O1+O2;
      O.name = ...
        sprintf('CNOT gate (%d qubits, control: q[%d], target: q[%d])',N,c,t);
    end

  end
  
  %% == Auxiliary Methods
  
  methods(Static)

    function disp(O)
      % Pretty-prints a qubit operator.
      %
      %   disp(O)
      %   disp(O,prec)
      %
      % arguments:
      %   O    - The operator, a Fock object.
      %   prec - Number of significant digits to print (optional, default is 4).
      %
      % returns:
      %   nothing

      if nargin<2; prec=4; end                                                  % Default number of sign. digits
      [~,N] = fockbasis(O).getNumSectors();                                     % Get number of qubits
      s = fock.HREF('ibmqx');                                                   % Hyperling to this class
      fprintf('  %s <strong>quantum gate operator</strong>\n\n',s);             % Print headline
      if ~strcmp(string(O.name),"")                                             % Non-empty name >>
        fprintf('    Name  : <strong>%s</strong>\n',O.name);                    %   Print name
      end                                                                       % <<
      fprintf('    Qubits: <strong>%d</strong>\n\n',N);                         % Print number of qubits
      fprintf(ibmqx.sprintop(O,'prec',prec));                                   % Print operator matrix
      fprintf('\n');                                                            % Print line break
    end

    function latex(O,prec)
      % Type-sets a qubit operator in LaTeX.
      %
      %   latex(O)
      %   latex(O,prec)
      %
      % arguments:
      %   O    - The operator, a Fock object.
      %   prec - Number of significant digits to print (optional, default is 4).
      %
      % returns:
      %   nothing
      
      if nargin<2; prec=4; end                                                  % Default number of sign. digits
      fprintf('%% Generated by FockBox'' ibmqx.latex method\n');                % LaTeX preamble
      fprintf(ibmqx.sprintop(O,'prec',prec,'mode','latex'));                    % Print operator matrix
      fprintf('\n');                                                            % Print line break
    end

    function s=sprintop(O,varargin)
      % Pretty-prints a quantum gate operator into a string.
      % 
      %   s=sprintop(O)
      %   s=sprintop(___,Name,Value)
      %
      % returns:
      %   s           - The string representation of the quantum gate
      %                 operator.
      %
      % arguments:
      %   O           - The operator, a Fock object.
      %   prec        - Number of significant digits to print (optional, default 
      %                 is 4).
      %   Name,Value  - Additional properties for printing.
      %   ___         - Any other argument list.
      %
      % additional properties:
      %   'prec'      - Number of significant digits to print. The default
      %                 value is 4.
      % 
      %   'mode'      - Printing mode, value can be one of the following
      %     'console' - For Matlab console (default)
      %     'latex'   - LaTeX code
      %     'plain'   - Plain ASCII

      % Initialize                                                              % -------------------------------------
      s = "";                                                                   % Initialize output string
      mode  = ibmqx.getModeOption(varargin);                                    % Get printing mode
      prec  = ibmqx.getPrecOption(varargin);                                    % Get no. of significant digits
      B     = fockbasis(O);                                                     % Get operator basis
      [~,N] = B.getNumSectors();                                                % Get number of qubits
      om    = ibmqx.opToMatrix(O);                                              % Get operator matrix

      % Preflight: Determine column widths                                      % -------------------------------------
      cw = zeros(size(om,1)+1,1);                                               % Initialize column with array
      for i=1:2^N                                                               % Loop over basis kets >>
        [~,l]   = ibmqx.sprintvec(ibmqx.int2binstr(i-1,N),'bra',mode);          %   Get printable length
        cw(1)   = max(cw(1),l);                                                 %   Aggregate with of header column
        cw(i+1) = l;                                                            %   Initialize with if header row
      end                                                                       % <<
      for i=1:size(om,1)                                                        % Loop over op. matrix rows >>
        for j=1:size(om,2)                                                      %   Loop over op. matrix columns >>
          [~,l] = ibmqx.sprintcmplx(om(i,j),prec,mode);                         %     Get printable len. of val. string
          cw(j+1) = max(cw(i+1),l);                                             %     Aggregate column width
        end                                                                     %   <<
      end                                                                       % <<

      % Pretty-print operator matrix                                            % -------------------------------------
      if strcmp(mode,'latex')                                                   % Latex mode >>
        sRH = '  '; sCD = ' & '; sRT = ' \\\\\n';                               %   Set table delimiters
        s = s + '%% Uses the kbordermatrix package!\n';                         %   ...
        s = s + '\\left(\\!\\kbordermatrix{\n';                                 %   ...
      else                                                                      % << All other modes >>
        sRH = '    '; sCD = ' '; sRT = '\n';                                    %   Set table delimiters
      end                                                                       % <<
      f = '%'+sprintf("%ds",cw(1)); s = s + sRH + sprintf(f,'');                % Output header of first row
      if strcmp(mode,'latex'); s = s + '    '; end                              % Extra LaTeX output
      for j=1:2^N                                                               % Loop over op. matrix columns >>
        s = s+sCD+ibmqx.sprintvec(ibmqx.int2binstr(j-1,N),'bra',mode,cw(j+1));  %   Output column header
      end                                                                       % <<
      s = s + sRT;                                                              % Output line break
      for i=1:size(om,1)                                                        % Loop over op. matrix rows >>
        s = s + sRH + ibmqx.sprintvec(ibmqx.int2binstr(i-1,N),'ket',mode,cw(1));%   Output row header
        if strcmp(mode,'latex'); s = s + '\\!\\!'; end                          %   Extra LaTeX output
        for j=1:size(om,2)                                                      %   Loop over op. matrix columns >>
          s = s + sCD + ibmqx.sprintcmplx(om(i,j),prec,mode,cw(j+1));           %     Output values
        end                                                                     %   <<
        if i<size(om,1); s = s + sRT; end                                       %   Output line break (except last row)
      end                                                                       % <<
      if strcmp(mode,'latex')                                                   % Latex mode >>
        s = s + '\n}\\!\\right)';                                               %   LaTeX postamble
      end                                                                       % <<
    end
    
  end  

  methods(Static,Access=private)

    function [s,len]=sprintcmplx(v,prec,mode,width)
      % Thrifty-prints a complex value to a formatted string (with markup).

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
        s = strrep(s,'i','\mathrm{i}');                                         %   Replace imaginary unit
        if ~strcmp(s,"0"); s = '\mathbf{'+s+'}'; end                            %   Print non-zero value bold
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
      
      if strcmp(type,'bra')                                                     % Print a bra vector >>
        skon = '<'; skoff = '|';                                                %   Set delimiters
        if strcmp(mode,'latex'); skon = '\langle '; end                         %   Set LaTeX left delimiter
      else                                                                      % << Print a ket vector >>
        skon = '|'; skoff = '>';                                                %   Set delimiters
        if strcmp(mode,'latex'); skoff = '\rangle'; end                         %   Set LaTeX right delimiter
      end                                                                       % <<
      if nargin<4; width=strlength(name); end                                   % Minimal print with
      f = "%"+sprintf('%ds',width);                                             % Make print format
      s = sprintf(f,string(skon)+string(name)+string(skoff));                   % Print ket
      len = strlength(s);                                                       % Printable length
      s = replace(s,'\','\\');                                                  % Escape backslash

    end

    function s=int2binstr(n,len)
      % Converts an non-negative integer to a binary number string.
      %
      %   s = int2binstr(n)
      %   s = int2binstr(n,len)
      %
      % arguments:
      %   n   - The integer.
      %   len - The length of the binary number.
      %
      % returns:
      %   s   - The binary number represented by a character array.

      if nargin<2; len=1; end
      s = "";
      while n>0
        b = rem(n,2);
        s = string(b) + s;
        n = floor(n/2);
      end
      while true
        if strlength(s)>=len; break; end
        s = "0"+s;
      end
      s = convertStringsToChars(s);
    end

    function mode=getModeOption(args)
      % Get pretty-printing mode from variable argument list.

      mode = 'console';                                                         % Default is Matlab console
      for i=1:2:length(args)                                                     % Loop over key-value pairs >>
        if strcmp(args{i},'mode')                                               %   Printing mode (normalize) >>
          if strcmp(args{i+1},'latex'); mode='latex'; end                       %     Latex
          if strcmp(args{i+1},'plain'); mode='plain'; end                       %     Plain ASCII
        end                                                                     %   <<
      end                                                                       % <<
    end

    function prec=getPrecOption(args)
      % Get number of significant digits for pretty-printing from variable 
      % argument list.
      
      prec = 4;                                                                 % Default is 4
      for i=1:2:length(args)                                                    % Loop over key-value pairs >>
        if strcmp(args{i},'prec')                                               %   No. sign. digits to print >>
          prec = args{i+1};                                                     %     Get value
        end                                                                     %   <<
      end                                                                       % <<
    end

  end

end

% EOF