function generateGrid

[f,v] = stlread('forma2.STL');
fv.faces = f;
fv.vertices = v;

figure
hold on;
patch(fv,'FaceColor','g','FaceAlpha',0.2)
axis equal

xmin = min(v(:,1));
xmax = max(v(:,1));
dx = xmax-xmin;
ymin = min(v(:,2));
ymax = max(v(:,2));
dy = ymax-ymin;
zmin = min(v(:,3));
zmax = max(v(:,3));
dz = zmax-zmin;

n = 5;
dl = dx/n;

n1 = n+1;
n2 = fix(dy/dl) + 1;
n3 = fix(dz/dl) + 1;

nc = (n1+1)*(n2+1)*(n3+1);
x = zeros(nc,1);
y = zeros(nc,1);
z = zeros(nc,1);

s = 1;
for i = 1:n1+1
    for j = 1:n2+1
        for k = 1:n3+1
            x(s) = xmin - dl/2 + (i-1)*dl;
            y(s) = ymin - dl/2 + (j-1)*dl;
            z(s) = zmin - dl/2 + (k-1)*dl;
            s = s + 1;
        end
    end
end

M = zeros(n1,n2,n3);
in = inpolyhedron(fv,[x,y,z]);

plot3(x,y,z,'.k')
plot3(x(in),y(in),z(in),'or')




