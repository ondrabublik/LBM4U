function [meshPath, simulPath, outputPath] = getPath(geometry, simulation)
LBMPath = getLBMPath;

if nargin == 0
    [geometry, simulation] = loadArgs;
end

meshPath = sprintf('%s/simulations/%s/%s/',LBMPath, geometry, 'mesh');
simulPath = sprintf('%s/simulations/%s/%s/',LBMPath, geometry, simulation);
outputPath = sprintf('%s/simulations/%s/%s/%s/',LBMPath, geometry, simulation, 'result');

if ~exist(meshPath, 'dir')
    error('Geometry %s does not exist.', geometry);
end

if ~exist(simulPath, 'dir')
    display('Simulation %s does not exist.');
    simulPath = [];
end

if ~exist(outputPath, 'dir')
    outputPath = [];
end
