function initBoundary
b = load('../mesh/blocks.txt');
ni = str2num(getValue('ni'));
nj = str2num(getValue('nj'));
nk = str2num(getValue('nk'));

nBi = max(b(:,1))+1;
nBj = max(b(:,2))+1;
mkdir('../mesh/boundary');
for bl = 1:length(b(:,1))
    iB = b(bl,1);
    jB = b(bl,2);
    kB = b(bl,3);
    indB = nBi * nBj * kB + nBi * jB + iB;
    fid = fopen(['../mesh/boundary/blockBoundary',num2str(indB),'.txt'],'w');
    
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                if(0)
                    fprintf(fid,'%d \n', 1);
                else
                    fprintf(fid,'%d \n', 0);
                end
            end
        end
    end
    fclose(fid);
    disp(['block: ',num2str(bl)]);
end
