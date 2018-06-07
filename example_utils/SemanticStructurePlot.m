classdef SemanticStructurePlot < handle
  % A graphviz semantic structure plot.
  %
  % authors:
  %   Markus Huber, BTU Cottbus-Senftenberg
  %
  % TODO:
  % -

  %% == Read-only properties ==

  properties(SetAccess=protected)

    % The axes of this plot.
    ax;

  end

  %% == Properties ==

  properties(SetAccess=private)

    % The value sets of the structure.
    vs;

    % The dependencies of the structure.
    ds;

    % The joint spaces of the structure.
    js;

  end

  %% == Public API ==

  methods

    function obj=SemanticStructurePlot()
      % Creates a semantic structure plot.
      %
      % returns:
      %   obj  - The semantic structure plot.
      %
      % example:
      %
      %   ssp = SemanticStructurePlot();

      % Create scene                                                            % -------------------------------------
      obj.ax = gca;                                                             % Remember axes
      axis equal tight off;                                                     % Parametrize plot axes

      obj.vs = containers.Map;                                                  % Create value sets map
      obj.ds = containers.Map;                                                  % Create dependencies map
      obj.js = containers.Map;                                                  % Create joint spaces map
    end

    function clear(obj)
      % Clears the plot.
      %
      %    obj.clear();
      %
      % arguments:
      %    obj - The semantic structure plot.

      img = findobj(obj.ax, 'Type', 'image');
      if ~isempty(img)
        delete(img);
      end

      obj.vs = containers.Map;                                                  % Create value sets map
      obj.ds = containers.Map;                                                  % Create dependencies map
      obj.js = containers.Map;                                                  % Create joint spaces map
    end

    function addValueSet(obj,name,vs,evs)
      % Adds a named value set to the structure.
      %
      %    obj.addValueSet(name,vs,evs);
      %
      % arguments:
      %    obj  - The semantic structure plot.
      %    name - The name for the value set.
      %    vs   - The fockobj representing all values (the sum of all values).
      %    evs  - The fockobj representing all seen values (again as sum).
      %
      % example:
      %
      %   obj.addValueSet('C', fockobj.bket('0C')+fockobj.bket('1C'), fockobj.bket('0C'));

      t.vs = vs;
      t.evs = evs;
      t.cn = 0;
      obj.vs(name) = t;
    end

    function addDependency(obj,ncon,cons)
      % Adds a dependency of a non-controllable on a controllable.
      %
      %    obj.addDependency(ncon,cons);
      %
      % arguments:
      %    obj  - The semantic structure plot.
      %    ncon - The name of the non-controllable.
      %    cons - The name(s) of the controllable(s) (a single name or a cell array of names).
      %
      % example:
      %
      %   obj.addDependency('C', {'X';'Y'});

      if (~obj.ds.isKey(ncon))
        obj.ds(ncon) = {};
      end

      t = obj.ds(ncon);
      t{end+1} = cons;
      obj.ds(ncon) = t;
    end

    function createGraph(obj,name)
      % Creates the graph as dot-file and png-file and displays it.
      %
      %    obj.createGraph(name);
      %
      % arguments:
      %    obj  - The semantic structure plot.
      %    name - The basename of the created files.
      %
      % example:
      %
      %   obj.createGraph('example');

      name = char(name);

      dot = fopen([name '.dot'], "w");
      fprintf(dot, 'graph {\n');

      wc = '';
      dsids = obj.ds.keys;
      for i=1:length(dsids)
        fprintf(dot, ['\n\tf_' dsids{i}...
                      '[label=<f<SUB><I>' dsids{i} '</I></SUB>>'...
                      ' texlbl="$f_{' dsids{i} '}$"'...
                      ' style="fspace/.try"];\n']);
        wc = [wc ' f_' dsids{i}];
        
        fprintf(dot, ['\tf_' dsids{i} ' -- { ' dsids{i}]);
        v = obj.vs(dsids{i});
        v.cn = 1;
        obj.vs(dsids{i}) = v;

        cons = obj.ds(dsids{i});
        for j=1:length(cons)
          t = cons{j};
          if (iscell(t))
            n = ['J_' strjoin(t,'')];
            fprintf(dot, [' ' n]);
            obj.js(n) = t;
          else
            fprintf(dot, [' ' t]);
            v = obj.vs(t);
            v.cn = 1;
            obj.vs(t) = v;
          end
        end
        fprintf(dot, ' };\n');
      end

      fprintf(dot, ['\n\t{ rank=same;' wc ' };\n']);

      jsids = obj.js.keys;
      for i=1:length(jsids)
        fprintf(dot, ['\n\t' jsids{i}...
                      '[label=<J<SUB><I>' strjoin(obj.js(jsids{i}),',') '</I></SUB>>'...
                      ' texlbl="$J_{' strjoin(obj.js(jsids{i}),',') '}$"'...
                      ' style="cspace/.try"];\n']);
        fprintf(dot, ['\t' jsids{i} ' -- {']);
        vsids = obj.js(jsids{i});
        for j=1:length(vsids)
          fprintf(dot, [' ' vsids{j}]);
          v = obj.vs(vsids{j});
          v.cn = 1;
          obj.vs(vsids{j}) = v;
        end
        fprintf(dot, ' };\n');
      end

      fprintf(dot, ['\n\tsubgraph Values {\n'...
                    '\t\tedge[dir=none, style=dotted]\n'...
                    '\t\tnode[shape=none]\n']);
      vns = '';
      vsids = obj.vs.keys;
      for i=1:length(vsids)
        vns = [vns ' ' vsids{i}...
               '[label=<<I>' vsids{i} '</I>>'...
               ' texlbl="$' vsids{i} '$"'...
               ' style="Vspace/.try"]'];
        v = obj.vs(vsids{i});
        if (~v.cn)
          wc = [wc ' ' vsids{i}];
        end
        rs = {};
        for j=1:size(v.vs.data,1)
          rs{end+1} = v.vs.data{j}(1:end-length(vsids{i}));
        end
        if (length(rs) < 4)
          rs = strjoin(rs, ',');
        else
          rs = [ rs{1} ',...,' rs{end} ];
        end
        ms = setdiff(v.vs.data(:,1),v.evs.data(:,1));
        if (~isempty(ms))
          if (length(ms) < 4)
            ms = ['\\{' strjoin(ms, ',') '}'];
          else
            ms = ['\\{' ms{1} ',...,' ms{end} '}'];
          end
        else
          ms = '';
        end
        fprintf(dot, ['\t\t' vsids{i} ' -- { ' vsids{i} 's'...
                      '[label=<{' rs '}' ms '>'...
                      ' texlbl="$\\{' texify(rs) '\\}' texify(ms) '$"'...
                      ' style="Vset/.try"] };\n']);
      end
      fprintf(dot, '\t}\n');

      fprintf(dot, ['\n\t{ rank=same;' vns ' }\n']);

      fprintf(dot, '\n\tW[texlbl="$\\mathcal{W}$" style="Wspace/.try"];\n');
      fprintf(dot, ['\tW -- {' wc ' };\n']);

      fprintf(dot, '}\n');
      fclose(dot);

      system(['dot -Tpng ' name '.dot > ' name '.png -q1']);
      img = imread([name '.png']);

      hold(obj.ax, 'on');
      imshow(img, 'Parent', obj.ax);
      hold(obj.ax, 'off');
      
      function tex=texify(str)
        tex = regexprep(str,'\\\','\\\\setminus');
        tex = regexprep(tex,'\.\.\.','\\\\dots');
        tex = regexprep(tex,'{','\\\\{');
        tex = regexprep(tex,'}','\\\\}');
      end
    end

  end

end

% EOF