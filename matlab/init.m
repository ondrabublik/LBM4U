function init
[meshPath, simulPath, outputPath] = getPath;
path = pwd;
cd(getLBMPath);
cd(simulPath);
try
    initBoundary;
catch
    'Boundary not defined!';
end

try
    initCondition;
catch
    'Init condition not defined!';
end
cd(path);