function generateGrid2

ni = 6;
nj = 6;
nk = 6;
deepCoarseLevel = 4;

[f,v] = stlread('ASMO.STL');

xmin = min(v(:,1));
xmax = max(v(:,1));
ymin = min(v(:,2));
zmin = min(v(:,3));
zmax = max(v(:,3));
dl = xmax-xmin;
v(:,1) = (v(:,1)-xmin)/dl;
v(:,2) = (v(:,2)-ymin)/dl;
v(:,3) = (v(:,3)-zmin)/dl - (zmax-zmin)/2/dl;
fv.faces = f;
fv.vertices = v;

figure
hold on;
patch(fv,'FaceColor','g','FaceAlpha',0.2)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

n1 = 20;
n2 = 7;
n3 = 8;
dl = 2/n1;

x0 = -0.5;
y0 = 1e-8;
z0 = -n3*dl/2;

nc = (n1+1)*(n2+1)*(n3+1);
x = zeros(nc,1);
y = zeros(nc,1);
z = zeros(nc,1);

s = 1;
for i = 1:n1+1
    for j = 1:n2+1
        for k = 1:n3+1
            x(s) = x0 + (i-1)*dl;
            y(s) = y0 + (j-1)*dl;
            z(s) = z0 + (k-1)*dl;
            s = s + 1;
        end
    end
end
in = inpolyhedron(fv,[x,y,z]);
ex = [0 1 1 0 0 1 1 0];
ey = [0 0 1 1 0 0 1 1];
ez = [0 0 0 0 1 1 1 1];

try
    which_dir = '/boundary';
    dinfo = dir(which_dir);
    dinfo([dinfo.isdir]) = [];   %skip directories
    filenames = fullfile(which_dir, {dinfo.name});
    delete( filenames{:} )
catch
    disp('Could not delete directory: /boundary')
end

disp('Buidling structure ...')
M = zeros(n1,n2,n3);
coarseLevel = ones(n1,n2,n3);
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            disp(['Block: ',num2str(n3*n2*(i-1) + n3*(j-1) + k-1)])
            M(i,j,k) = 1;
            b = 0;
            for p = 1:8
                ind = (n3+1)*(n2+1)*(i-1 + ex(p)) + (n3+1)*(j-1 + ey(p)) + k + ez(p);
                if(in(ind) == 1)
                    b = 1;
                end
            end
            if(b == 1)
                coarseLevel(i,j,k) = deepCoarseLevel;
                nib = ni*coarseLevel(i,j,k);
                njb = nj*coarseLevel(i,j,k);
                nkb = nk*coarseLevel(i,j,k);
                
                nct = nib*njb*nkb;
                xt = zeros(nct,1);
                yt = zeros(nct,1);
                zt = zeros(nct,1);

                r = (n3+1)*(n2+1)*(i-1) + (n3+1)*(j-1) + k;
                s = 1;
                for it = 1:nib
                    for jt = 1:njb
                        for kt = 1:nkb
                            xt(s) = x(r) + (it-1)/(nib-1)*dl;
                            yt(s) = y(r) + (jt-1)/(njb-1)*dl;
                            zt(s) = z(r) + (kt-1)/(nkb-1)*dl;
                            s = s + 1;
                        end
                    end
                end
                int = inpolyhedron(fv,[xt,yt,zt]);
                plot3(xt(int),yt(int),zt(int),'.g')
                
                bl = n3*n2*(i-1) + n3*(j-1) + k;
                fid = fopen(['../../boundary/blockBoundary',num2str(bl-1),'.txt'],'w');
                for it = 1:nib
                    for jt = 1:njb
                        for kt = 1:nkb
                            b = 1;
                            ind = njb*nkb*(it-1) + nkb*(jt-1) + kt;
                            if(int(ind) == 1)
                                b = 0;
                            end
                            fprintf(fid,'%i\n',b);
                        end
                    end
                end
                fclose(fid);
            end
        end
    end
end
% plot3(x,y,z,'.k')
% plot3(x(in),y(in),z(in),'or')

ex = [-1 1 0 0 0 0 -1 1 -1 1 0 0 0 0 -1 1 -1 1 -1 -1 -1 -1 1 1 1 1];
ey = [0 0 -1 1 0 0 -1 -1 1 1 -1 1 -1 1 0 0 0 0 -1 -1 1 1 -1 -1 1 1];
ez = [0 0 0 0 -1 1 0 0 0 0 -1 -1 1 1 -1 -1 1 1 -1 1 -1 1 -1 1 -1 1];
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            for p = 1:26
                try
                    if(coarseLevel(i,j,k) == 1 && coarseLevel(i+ex(p),j+ey(p),k+ez(p)) == 4)
                        coarseLevel(i,j,k) = 2;
                    end
                catch
                end
            end
        end
    end
end

generateStructure(M,coarseLevel,1);
setValue('../../parameters.txt','ni',num2str(ni));
setValue('../../parameters.txt','nj',num2str(nj));
setValue('../../parameters.txt','nk',num2str(nk));




