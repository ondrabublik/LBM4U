function path = getLBMPath
    fid = fopen('LBMPath.txt','r');
    path = fgetl(fid);
    fclose(fid);