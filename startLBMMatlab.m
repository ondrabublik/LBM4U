
function startLBMMatlab
path = pwd;
addpath([path,'/matlab']);
fid = fopen('matlab/LBMPath.txt','w');
fprintf(fid,'%s\n',path);
fclose(fid);
