function coarseHexa = genCoarse

M = ones(5,5,5);
ne = sum(sum(sum(M))); % number of filled blocks
coarseHexa{ne} = '';
s = 1;
for i = 1:size(M,1)
    for j = 1:size(M,2)
        for k = 1:size(M,3)
            if(M(i,j,k) == 1)
                for i1 = 1:2
                    for j1 = 1:2
                        for k1 = 1:2
                            bin = '';
                            for i2 = 1:2
                                for j2 = 1:2
                                    for k2 = 1:2
                                        if(zahusti(i-1 + (i1-1)/2 + (i2-1)/4, j-1 + (j1-1)/2 + (j2-1)/4, k-1 + (k1-1)/2 + (k2-1)/4))
                                            bin = [bin,'1 '];
                                        else
                                            bin = [bin,'0 '];
                                        end
                                    end
                                end
                            end
                            coarseHexa{s} = [coarseHexa{s},' ',toHexa(bin)];
                        end
                    end
                end
                s = s + 1;
            end
        end
    end
end



%================================================================
% figure
% hold on
% for i = 1:length(M(:,1,1))
%     for j = 1:length(M(1,:,1))
%         for k = 1:length(M(1,1,:))
%             if(M(i,j,k) == 1)
%                 drawBox(i,j,k);
%             end
%         end
%     end
% end
% axis equal
% 
% name = inputdlg('Geometry name: ');
% name = ['../../simulations/',name{1}];
% if(isdir(name))
%     answer = questdlg('Geometry exists. Rewrite?');
%     if(~answer)
%         return;
%     end  
% else
%     mkdir(name);
% end
% name = [name,'/mesh'];
% 
% mkdir(name);
% fid = fopen([name,'/blocks.txt'],'w');
% for i = 1:length(M(:,1,1))
%     for j = 1:length(M(1,:,1))
%         for k = 1:length(M(1,1,:))
%             if(M(i,j,k) == 1)
%                 fprintf(fid,'%d %d %d\n',i-1,j-1,k-1);
%             end
%         end
%     end
% end
% fclose(fid);
% close all;
% 
% 
% function drawBox(i,j,k)
% x = [0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 0];
% y = [0 0 1 1 0 0 0 0 0 1 1 1 1 1 1 0];
% z = [0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1];
% plot3(x + i,z+j,y+k,'color','k','linewidth',2)
% xlabel('X','fontsize',16);
% ylabel('Y','fontsize',16);
% zlabel('Z','fontsize',16);


function out = toHexa(in)
    terms = strsplit(in);
    count = 0;
    for i = 1:8
        count = count + str2num(terms{i})*2^(8-i);
    end
    out = dec2hex(count);
    if(length(out) == 1)
        out = ['0',out];
    end

function out = zahusti(i,j,k)
r = sqrt((i-2.5)^2 + (j-2.5)^2 + (k-2.5)^2);
out = 0;
if(r < 1)
    out = 1;
end





