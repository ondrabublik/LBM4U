function generateBoundary
global xmax ymax zmax
gridPos = load('../../blockGlobalGridPosition.txt');
levels = load('../../blocksLevel.txt');
ni = str2num(getValue('../../parameters.txt','ni'));
nj = str2num(getValue('../../parameters.txt','nj'));
nk = str2num(getValue('../../parameters.txt','nk'));
dx = 1/max(gridPos(:,1)+1);
xmax = max(gridPos(:,1))*dx;
ymax = max(gridPos(:,2))*dx;
zmax = max(gridPos(:,3))*dx;

nBlock = size(gridPos,1);

for bl = 1:nBlock
    disp(['block: ', num2str(bl)]);
    fid = fopen(['../../boundary/blockBoundary',num2str(bl-1),'.txt'],'w');
    nib = ni*levels(bl);
    njb = nj*levels(bl);
    nkb = nk*levels(bl);
    for i = 1:nib
        for j = 1:njb
            for k = 1:nkb
                b = 1;
                x = (gridPos(bl,1)-1)*dx + (i-1.5)/(nib-1)*dx;
                y = (gridPos(bl,2)-1)*dx + (j-1.5)/(njb-1)*dx;
                z = (gridPos(bl,3)-1)*dx + (k-1.5)/(nkb-1)*dx;
                if(inShape(x,y,z))
                    b = 0;
                end
                fprintf(fid,'%i\n',b);
            end
        end
    end
    fclose(fid);
end

function in = inShape(x,y,z)
global xmax ymax zmax
in = 0;
r = sqrt((x-xmax/2)^2 + (y-ymax/2)^2 + (z-zmax/2)^2);
if(r < xmax/15)
    in = 1;
end

