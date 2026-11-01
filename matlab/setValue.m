function setValue(propName,propVal)

[meshPath, simulPath, outputPath] = getPath;
fileName = [simulPath,'parameters.txt'];

% Read txt into cell A
fid = fopen(fileName,'r');
tline = fgetl(fid);
A{1} = [];
i = 1;
while ischar(tline)
    A{i} = tline;
    str = strsplit(tline,'=');
    str = delWhitespace(str{1});
    if(strcmp(str,propName))
        A{i} = [propName,' = ', propVal];
    end
    tline = fgetl(fid);
    i = i + 1;
end
fclose(fid);

fid = fopen(fileName, 'w');
for i = 1:length(A)
    fprintf(fid,'%s\n', A{i});
end
fclose(fid);

function str = delWhitespace(str)
    is = 1;
    ik = length(str);
    
    while(str(is) == ' ')
        is = is + 1;
    end
    
    while(str(ik) == ' ')
        ik = ik - 1;
    end
    str = str(is:ik);