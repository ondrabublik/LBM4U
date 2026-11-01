function test

M = zeros(5,3,3);
M(1:5,1:3,1:3) = 1;

coarseLevel = ones(5,3,3);
for i = 1:5
    for j = 1:3
        for k = 1:3
            r = sqrt((i-3)^2 + (j-2)^2 + (k-2)^2);
            if(r < 1)
               coarseLevel(i,j,k) = 4;
            elseif(r < 2)
               coarseLevel(i,j,k) = 2;
            end
        end
    end
end

M = ones(1,1,1);
coarseLevel = ones(size(M));

generateStructure(M,coarseLevel,1);
generateBoundary;