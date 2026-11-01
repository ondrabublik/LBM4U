function generateStructure(M,coarseLevel,periodX)

nBlock = sum(sum(sum(M)));
[n1,n2,n3] = size(M);
blockStruct = -1*ones(nBlock,26);
blockIndex = ones(nBlock,3);

ex = [-1 1 0 0 0 0 -1 1 -1 1 0 0 0 0 -1 1 -1 1 -1 -1 -1 -1 1 1 1 1];
ey = [0 0 -1 1 0 0 -1 -1 1 1 -1 1 -1 1 0 0 0 0 -1 -1 1 1 -1 -1 1 1];
ez = [0 0 0 0 -1 1 0 0 0 0 -1 -1 1 1 -1 -1 1 1 -1 1 -1 1 -1 1 -1 1];

blockName = zeros(n1,n2,n3);
s = 0;
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            if(M(i,j,k) == 1)
                blockName(i,j,k) = s;
                s = s + 1;
            end
        end
    end
end

s = 1;
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            if(M(i,j,k) == 1)
                blockIndex(s,:) = [i,j,k];
                for p = 1:26
                    ip = i+ex(p);
                    if(periodX)
                        if(ip > n1)
                            ip = 1;
                        end
                        if(ip < 1)
                            ip = n1;
                        end
                    end
                    try
                        if(M(ip,j+ey(p),k+ez(p)) == 1)
                            blockStruct(s,p) = blockName(ip,j+ey(p),k+ez(p));
                        end
                    catch
                    end
                end
                s = s + 1;
            end
        end
    end
end

fid = fopen('../../blocks.txt','w');
for i = 1:size(blockStruct,1)
    for j = 1:size(blockStruct,2)
        fprintf(fid,'%d ',blockStruct(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('../../blockGlobalGridPosition.txt','w');
for i = 1:size(blockIndex,1)
    for j = 1:size(blockIndex,2)
        fprintf(fid,'%d ',blockIndex(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('../../blocksLevel.txt','w');
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            if(M(i,j,k) == 1)
                fprintf(fid,'%d ',coarseLevel(i,j,k));
                fprintf(fid,'\n');
            end
        end
    end
end
fclose(fid);

figure
hold on
for i = 1:nBlock
    drawBox(blockIndex(i,:),num2str(i-1));
end
axis equal



function drawBox(ind,str)
x = [0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 0];
y = [0 0 1 1 0 0 0 0 0 1 1 1 1 1 1 0];
z = [0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1];
plot3(x+ind(1),z+ind(3),y+ind(2),'color','k','linewidth',2)
text(ind(1),ind(3),ind(2),str,'fontsize',16,'color','r')
xlabel('X','fontsize',16);
ylabel('Z','fontsize',16);
zlabel('Y','fontsize',16);








