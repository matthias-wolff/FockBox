classdef SemanticStructurePlot < handle
  % A semantic structure plot.
  %
  % authors:
  %   Markus Huber, BTU Cottbus-Senftenberg
  %
  % TODO:
  % - add value sets to plot
  % - use adjacency matrix for creation of dot-file

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

    % The connection vector of W.
    cvs;

    % The current plot objects.
    po;

    % The adjacency matrix of the graph.
    aj;
  
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
      set(gcf,'color','w');                                                     % Set figure plot background to white

      obj.clear();
    end

    function clear(obj)
      % Clears the plot.
      %
      %    obj.clear();
      %
      % arguments:
      %    obj - The semantic structure plot.

      obj.vs = containers.Map;                                                  % Create value sets map
      obj.ds = containers.Map;                                                  % Create dependencies map
      obj.js = containers.Map;                                                  % Create joint spaces map
      obj.cvs = zeros(0);                                                       % Create connection vector of W

      obj.preview();
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
      t.in = length(obj.vs)+1;
      t.n = ['$' name '$'];
      obj.vs(name) = t;

      t = zeros(1, length(obj.cvs)+1);
      t(1, 1:end-1) = obj.cvs;
      t(1, end) = 1;
      obj.cvs = t;
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
      %   obj.addDependency('C', {'X','Y'});

      if (~obj.ds.isKey(ncon))
        t.in = length(obj.ds)+1;
        t.n = ['$f_{' ncon '}$'];
        t.cs = {};
        t.cvs = zeros(1, length(obj.vs));
        t.cvs(1, obj.vs(ncon).in) = 1;
        obj.cvs(1, obj.vs(ncon).in) = 0;
        t.cjs = zeros(0);
      else
        t = obj.ds(ncon);
      end

      if iscell(cons)
        if length(cons) > 1
          n = ['J_' strjoin(cons, '')];
          if obj.js.isKey(n)
            jid = obj.js(n).in;
            if jid < length(t.cjs)
              t.cjs(1, jid) = 1;
            else
              j = zeros(1, jid);
              j(1, 1:length(t.cjs)) = t.cjs;
              j(1, jid) = 1;
              t.cjs = j;
            end
          else
            s.in = length(obj.js)+1;
            s.n = ['$J_{' strjoin(cons, ',') '}$'];
            s.cvs = zeros(1, length(obj.vs));
            for i=1:length(cons)
              s.cvs(1, obj.vs(cons{i}).in) = 1;
              obj.cvs(1, obj.vs(cons{i}).in) = 0;
            end
            obj.js(n) = s;

            j = zeros(1, length(obj.js));
            j(1, 1:length(t.cjs)) = t.cjs;
            j(1, s.in) = 1;
            t.cjs = j;
          end
        else
          t.cvs(1, obj.vs(cons{1}).in) = 1;
          obj.cvs(1, obj.vs(cons{1}).in) = 0;
        end
      else
        t.cvs(1, obj.vs(cons).in) = 1;
        obj.cvs(1, obj.vs(cons).in) = 0;
      end

      t.cs{end+1} = cons;
      obj.ds(ncon) = t;
    end

    function createGraph(obj,name)
      % Creates the graph as dot-file.
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
      fprintf(dot, 'graph G {\n');

      wc = '';
      dsids = obj.ds.keys;
      for i=1:length(dsids)
        fprintf(dot, ['\n\tf_' dsids{i}...
                      '[label=<f<SUB><I>' dsids{i} '</I></SUB>>'...
                      ' texlbl="$f_{' dsids{i} '}$"'...
                      ' style="fspace/.try"];\n']);
        wc = [wc ' f_' dsids{i}];
        
        fprintf(dot, ['\tf_' dsids{i} ' -- { ' dsids{i}]);

        d = obj.ds(dsids{i});
        cons = d.cs;
        for j=1:length(cons)
          t = cons{j};
          if (iscell(t))
            n = ['J_' strjoin(t,'')];
            fprintf(dot, [' ' n]);
            obj.js(n) = t;
          else
            fprintf(dot, [' ' t]);
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
                      ' style="jspace/.try"];\n']);
        fprintf(dot, ['\t' jsids{i} ' -- {']);
        vsids = obj.js(jsids{i});
        for j=1:length(vsids)
          fprintf(dot, [' ' vsids{j}]);
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
        if obj.aj(1,1+v.in)
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
            es = {};
            for j=1:length(ms)
              es{end+1} = ms{j}(1:end-length(vsids{i}));
            end
            ms = ['\\{' strjoin(es, ',') '}'];
          else
            a = ms{1}(1:end-length(vsids{i}));
            e = ms{end}(1:end-length(vsids{i}));
            ms = ['\\{' a ',...,' e '}'];
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

      function tex=texify(str)
        tex = regexprep(str, '\\\', '\\\\setminus');
        tex = regexprep(tex, '\.\.\.', '\\\\dots');
        tex = regexprep(tex, '{', '\\\\{');
        tex = regexprep(tex, '}', '\\\\}');
      end
    end

    function computeAdjacencyMatrix(obj)
      ovs = 2;
      lvs = length(obj.vs);
      ods = 2+lvs;
      lds = length(obj.ds);
      ojs = 2+lvs+lds;
      ljs = length(obj.js);

      t = zeros(1+lvs+lds+ljs);

      t(1,ovs:(ovs-1)+lvs) = obj.cvs;
      t(ovs:(ovs-1)+lvs,1) = obj.cvs';
      t(1,ods:(ods-1)+lds) = ones(1,lds);
      t(ods:(ods-1)+lds,1) = ones(lds,1);

      dsids = obj.ds.keys;
      for i=1:length(dsids)
        d = obj.ds(dsids{i});
        lvs = length(d.cvs);
        t((ods-1)+d.in,ovs:(ovs-1)+lvs) = d.cvs;
        t(ovs:(ovs-1)+lvs,(ods-1)+d.in) = d.cvs';
        ljs = length(d.cjs);
        t((ods-1)+d.in,ojs:(ojs-1)+ljs) = d.cjs;
        t(ojs:(ojs-1)+ljs,(ods-1)+d.in) = d.cjs';
      end

      jsids = obj.js.keys;
      for i=1:length(jsids)
        j = obj.js(jsids{i});
        lvs = length(j.cvs);
        t((ojs-1)+j.in,ovs:(ovs-1)+lvs) = j.cvs;
        t(ovs:(ovs-1)+lvs,(ojs-1)+j.in) = j.cvs';
      end

      obj.aj = t;
    end
    
    function preview(obj)
      if isstruct(obj.po)
        delete(obj.po.plot);
        delete(obj.po.nodes);
      end

      po.plot = gobjects(1,2);
      po.nodes = gobjects(1, 1+length(obj.vs)+length(obj.ds)+length(obj.js));

      obj.computeAdjacencyMatrix();

      hold(obj.ax, 'on');
      if length(obj.vs)
        po.plot(1) = plot(obj.ax, graph(obj.aj), 'k-o', 'Layout', 'layered',...
                          'NodeLabel', '', 'MarkerSize', 30, 'Linewidth', 2,...
                          'Sources', 1, 'Sinks', [2:length(obj.vs)+1]);
        po.plot(2) = plot(obj.ax, graph(obj.aj), 'ko', 'Layout', 'layered',...
                          'NodeLabel', '', 'MarkerSize', 28, 'Linewidth', 2,...
                          'NodeColor', 'w',...
                          'Sources', 1, 'Sinks', [2:length(obj.vs)+1]);
      else
        po.plot(1) = plot(obj.ax, graph(obj.aj), 'k-o', 'Layout', 'layered',...
                          'NodeLabel', '', 'MarkerSize', 30, 'Linewidth', 2,...
                          'Sources', 1);
        po.plot(2) = plot(obj.ax, graph(obj.aj), 'ko', 'Layout', 'layered',...
                          'NodeLabel', '', 'MarkerSize', 28, 'Linewidth', 2,...
                          'NodeColor', 'w',...
                          'Sources', 1);
      end
      hold(obj.ax, 'off');

      o = 0;
      id = o+1;
      lb = '$\mathcal{W}$';
      x = po.plot(1).XData(id);
      y = po.plot(1).YData(id);
      po.nodes(id) = text(obj.ax, x, y, lb,...
                          'Interpreter', 'latex',...
                          'HorizontalAlignment', 'center');
      o = o+1;

      vsids = obj.vs.keys;
      for i=1:length(vsids)
        id = o+obj.vs(vsids{i}).in;
        lb = obj.vs(vsids{i}).n;
        x = po.plot(1).XData(id);
        y = po.plot(1).YData(id);
        po.nodes(id) = text(obj.ax, x, y, lb,...
                            'Interpreter', 'latex',...
                            'HorizontalAlignment', 'center');
      end
      o = o+i;

      dsids = obj.ds.keys;
      for i=1:length(dsids)
        id = o+obj.ds(dsids{i}).in;
        lb = obj.ds(dsids{i}).n;
        x = po.plot(1).XData(id);
        y = po.plot(1).YData(id);
        po.nodes(id) = text(obj.ax, x, y, lb,...
                            'Interpreter', 'latex',...
                            'HorizontalAlignment', 'center');
      end
      o = o+i;


      jsids = obj.js.keys;
      for i=1:length(jsids)
        id = o+obj.js(jsids{i}).in;
        lb = obj.js(jsids{i}).n;
        x = po.plot(1).XData(id);
        y = po.plot(1).YData(id);
        po.nodes(id) = text(obj.ax, x, y, lb,...
                            'Interpreter', 'latex',...
                            'HorizontalAlignment', 'center');
      end

      obj.po = po;

      % - Redraw                                                                % - - - - - - - - - - - - - - - - - - -
      drawnow;                                                                  % Update plots  
    end

  end

end

% EOF