function addFockBoxToPath()
  addpath(strrep(mfilename('fullpath'),mfilename,'..'));                        % Add FockBox folder to path
  addpath(strrep(mfilename('fullpath'),mfilename,'../example_utils'));          % Add example utilities folder to path
end
