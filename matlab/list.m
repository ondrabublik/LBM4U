function [geomID,simID] = list(ID)

geomID = '';
simID = '';
if(nargin == 0)
    ID = -1;
end

LBMPath = getLBMPath;
geomList = dir([LBMPath,'/simulations']);
s = 1;
for i = 1 : length(geomList)
    if geomList(i).isdir && ~strcmp(geomList(i).name, '.') && ~strcmp(geomList(i).name , '..')
        simulList = dir(strcat(LBMPath,'/simulations/', geomList(i).name));
        for j = 1 : length(simulList)
            if simulList(j).isdir && ~strcmp(simulList(j).name, '.') && ~strcmp(simulList(j).name , '..') && ~strcmp(simulList(j).name , 'mesh')
                if(ID == -1)
                    fprintf(1, '%d %s %s \n', s, geomList(i).name, simulList(j).name);
                end
                if(s == ID)
                    geomID = geomList(i).name;
                    simID = simulList(j).name;
                end
                s = s + 1;
            end
        end
    end
end
