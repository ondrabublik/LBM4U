function run
clc
path = pwd;
cd(getLBMPath);
try
    %copyfile('dist\Release\MinGW_1-Windows\clbmnet.exe','clbm.exe')
    copyfile('vs\x64\Release\cLBMnet.exe','clbm.exe')
catch
    disp('exe file not found!')
end
tic
% system('cLBMnet.exe')
system('clbm.exe');
% system('start clbm.exe');
toc
cd(path)
