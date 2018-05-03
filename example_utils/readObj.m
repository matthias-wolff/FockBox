function obj = readObj(objFname,arg2)
  % Reads an wavefront object file.
  % This function parses wavefront object data. It reads the mesh vertices,
  % texture coordinates, normal coordinates and face definitions(grouped by
  % number of vertices) in a .obj file.
  %
  %   obj = readObj(objFname)
  %   obj = readObj(objFname,texFname)
  %   obj = readObj(objFname,[R G B])
  %
  % arguments:
  %    objFname - Wavefront object file full path
  %    texFname - Texture file full path
  %    [R G B]  - Object color
  %
  % returns:
  %    obj.v    - Mesh vertices
  %    obj.vt   - Texture coordinates
  %    obj.vn   - Normal coordinates
  %    obj.f    - Face definition assuming faces are made of of 3 vertices
  %    obj.g    - Group names and index of first face in group
  %    obj.fvcd - Face vertex color data
  %
  % author:
  %    Bernard Abayowa, Tec^Edge, 11/8/07
  %    Matthias Wolff, BTU Cottbus-Senftenberg (corrections and additions)
  %    - support of negative vertex indices in face definitions
  %    - suppert of group tier
  %
  % references:
  % [1] https://de.mathworks.com/matlabcentral/fileexchange/18957-readobj
  %     visited Mar. 31, 2018
  % [2] https://de.mathworks.com/matlabcentral/fileexchange/20307-display-obj
  %     visited Mar. 31, 2018

  % set up field types
  v = []; vt = []; vn = []; f.v = []; f.vt = []; f.vn = []; g = {};

  fid = fopen(objFname);

  % parse .obj file
  lctr = 0;
  vctr = 1;
  vnctr = 1;
  vtctr = 1;
  fctr = 1;
  gctr = 1;
  while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file
    ln = sscanf(tline,'%s',1); % line type
    lctr = lctr+1;
    %disp(ln)
    switch ln
      case 'g'   % Group name
        g{gctr} = {tline(3:end) fctr};
        gctr = gctr+1;
      case 'v'   % mesh vertexs
        v = [v; sscanf(tline(2:end),'%f')'];
        vctr = vctr+1;
      case 'vt'  % texture coordinate
        vt = [vt; sscanf(tline(3:end),'%f')'];
        vtctr = vtctr+1;
      case 'vn'  % normal coordinate
        vn = [vn; sscanf(tline(3:end),'%f')'];
        vnctr = vnctr+1;
      case 'f'   % face definition
        fv = []; fvt = []; fvn = [];
        str = textscan(tline(2:end),'%s'); str = str{1};

        nf = length(findstr(str{1},'/')); % number of fields with this face vertices

        [tok str] = strtok(str,'//');     % vertex only
        for k = 1:length(tok) fv = [fv str2num(tok{k})]; end

        if (nf > 0)
          [tok str] = strtok(str,'//');   % add texture coordinates
          for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
        end
        if (nf > 1)
          [tok str] = strtok(str,'//');   % add normal coordinates
          for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
        end
        % Handle negative vertex ids
        for i=1:length(fv);  if fv (i)<0; fv (i)=vctr + fv(i); end; end
        for i=1:length(fvt); if fvt(i)<0; fvt(i)=vtctr+fvt(i); end; end
        for i=1:length(fvn); if fvn(i)<0; fvn(i)=vnctr+fvn(i); end; end

        % Adjust array sizes for concatenation
        while size(fv ,2)<size(f.v ,2); fv  = [fv  NaN]; end
        while size(fvt,2)<size(f.vt,2); fvt = [fvt NaN]; end
        while size(fvn,2)<size(f.vt,2); fvn = [fvn NaN]; end

        while size(f.v,2)<size(fv,2); f.v = [f.v repmat(NaN,size(f.v,1),1)]; end
        while size(f.vt,2)<size(fvt,2); f.vt = [f.vt repmat(NaN,size(f.vt,1),1)]; end
        while size(f.vn,2)<size(fvn,2); f.vn = [f.vn repmat(NaN,size(f.vn,1),1)]; end

        f.v = [f.v; fv];
        f.vt = [f.vt; fvt];
        f.vn = [f.vn; fvn];
        fctr = fctr+1;
    end
  end
  fclose(fid);

  % Set up matlab object
  obj.v = v; obj.vt = vt; obj.vn = vn; obj.f = f; obj.g = g;
  
  % Make face vertex color data
  if nargin>1
    if ischar(arg2)||isstring(arg2)
      % arg2 is a texture file name
      texture_img = imread(arg2);
      [sy,sx,sz] = size(texture_img);
      texture_img = reshape(texture_img,sy*sx,sz);
      
      % Convert greyscale to RGB
      if sz == 1
        texture_img = repmat(texture_img,1,3);
      end
      
      [~,fv_idx] = unique(obj.f.v);
      texture_idx = obj.f.vt(fv_idx);
      texture_idx(isnan(texture_idx)) = [];
      x = abs(round(obj.vt(:,2)*(sx-1)))+1;
      y = abs(round(obj.vt(:,1)*(sy-1)))+1;
      xy = sub2ind([sy sx],y,x);
      texture_pts = xy(texture_idx);
      obj.fvcd = double(texture_img(texture_pts,:))/255;
    else
      % arg2 is an RGB triplet
      obj.fvcd = repmat(arg2,size(obj.v,1),1);
    end
  else
    % arg2 is empty
    obj.fvcd = repmat([0.5 0.5 0.5],size(obj.v,1),1);
  end

end

% EOF