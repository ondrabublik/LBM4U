%param  Open text file with parameters for the problem.

function param

[geometry, simulation] = loadArgs;
LBMPath = getLBMPath;
actPath = pwd;
cd(LBMPath);
path = sprintf('%s/simulations/%s/', LBMPath, geometry);
if ~exist(path, 'dir')
    error('Geometry %s does not exist.', geometry);
end
path = sprintf('%s%s/', path, simulation);
if ~exist(path, 'dir')
    answer = input(sprintf('Simulation %s does not exist, do you want to creat it [Y/n]? ', ...
             simulation), 's');
    
    switch lower(answer)
        case {'y', 'yes', ''}
            mkdir(path);
            copyfile('parameters.txt',[path,'parameters.txt']);
            copyfile('matlab/initCondition.m',[path,'initCondition.m']);
            copyfile('matlab/initBoundary.m',[path,'initBoundary.m']);
        otherwise
            return;
    end
end

path = strcat(path, 'parameters.txt');
edit(path);
cd(actPath);

jumpToSim;