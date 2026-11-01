function genStruct(save)

if(nargin == 0)
    save = 0;
end

M = ones(12,5,1);
M(4,2:end,1) = 0;
M(6,2:end,1) = 0;
M(9,2:end,1) = 0;
M(1:3,4:end) = 0;
M(10:12,4:end) = 0;


% M(1:11,1:3,1) = 1;
% M(1:5,5:9,1) = 1;
% M(1,4,1) = 1;
% M(1:3,1,1) = 0;
% M(6:8,1,1) = 0;
% M([1,3,6,8],2,1) = 0;


%================================================================
figure
hold on
for i = 1:length(M(:,1,1))
    for j = 1:length(M(1,:,1))
        for k = 1:length(M(1,1,:))
            if(M(i,j,k) == 1)
                drawBox(i,j,k);
            end
        end
    end
end
axis equal

if(save)
    name = inputdlg('Geometry name: ');
    name = ['../../simulations/',name{1}];
    if(isdir(name))
        answer = questdlg('Geometry exists. Rewrite?');
        if(~answer)
            return;
        end  
    else
        mkdir(name);
    end
    name = [name,'/mesh'];

    mkdir(name);
    fid = fopen([name,'/blocks.txt'],'w');
    for i = 1:length(M(:,1,1))
        for j = 1:length(M(1,:,1))
            for k = 1:length(M(1,1,:))
                if(M(i,j,k) == 1)
                    fprintf(fid,'%d %d %d\n',i-1,j-1,k-1);
                end
            end
        end
    end
    fclose(fid);
end


function drawBox(i,j,k)
x = [0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 0];
y = [0 0 1 1 0 0 0 0 0 1 1 1 1 1 1 0];
z = [0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1];
plot3(x + i,z+j,y+k,'color','k','linewidth',2)
xlabel('X','fontsize',16);
ylabel('Y','fontsize',16);
zlabel('Z','fontsize',16);


