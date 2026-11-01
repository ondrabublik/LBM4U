function initCondition
b = load('../mesh/blocks.txt');
ni = str2num(getValue('ni'));
nj = str2num(getValue('nj'));
nk = str2num(getValue('nk'));

nBi = max(b(:,1))+1;
nBj = max(b(:,2))+1;

mkdir('init');
for bl = 1:length(b(:,1))
    iB = b(bl,1);
    jB = b(bl,2);
    kB = b(bl,3);
    indB = nBi * nBj * kB + nBi * jB + iB;
    fid  = fopen(['init/blockInit',num2str(indB),'.txt'],'w');
    fid2 = fopen(['init/concentration',num2str(indB),'.txt'],'w');
    
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                % RHO,U,V,W,FLAG
%                 if (((i-20)^2+(j-20)^2+(k-20)^2) <= 10^2)
                    fprintf(fid,'%6.4f %6.4f %6.4f %6.4f %d \n', 1,0,0,0,1);
%                 else
%                     fprintf(fid,'%6.4f %6.4f %6.4f %6.4f %d \n', 1,0,0,0,2);
%                 end
                % CONCENTRATION
                if (((i-20)^2+(j-20)^2+(k-20)^2) <= 10^2)
                    fprintf(fid2,'%6.4f \n', 0);
                else
                    fprintf(fid2,'%6.4f \n', 0);
                end
            end
        end
    end
    fclose(fid);
    fclose(fid2);
    disp(['block: ',num2str(bl)]);
end

