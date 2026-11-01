function paraShow(varargin)
[~, simulPath, outputPath] = getPath;
if nargin == 0
    state = '';
    file = 'output_..us.vtk';
end

if nargin == 1
    state = '';
    file = varargin{1};
end

if nargin == 2
    state = ['--state=',simulPath,varargin{1}];
    file = varargin{2};
end
['paraview ', state, ' ', outputPath, file]
system(['paraview.exe ', state, ' ', outputPath, file]);