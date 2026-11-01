function generateGrid3

ni = 15;
nj = 15;
nk = 15;

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

n1 = 5;
n2 = 2;
n3 = 4;
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

mkdir('boundary');

disp('Buidling structure ...')
fid0 = fopen('blocks.txt','w');
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            b = 0;
            for p = 1:8
                ind = (n3+1)*(n2+1)*(i-1 + ex(p)) + (n3+1)*(j-1 + ey(p)) + k + ez(p);
                if(in(ind) == 0)
                    b = 1;
                end
            end
            if(b == 1)
                disp(['Block: ',num2str(n3*n2*(i-1) + n3*(j-1) + k-1)])
                fprintf(fid0,'%d %d %d\n',i-1,j-1,k-1);
                nib = ni;
                njb = nj;
                nkb = nk;
                
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
                plot3(xt(int==0),yt(int==0),zt(int==0),'.r')
                
                bl = n3*n2*(i-1) + n3*(j-1) + k;
                fid = fopen(['boundary/blockBoundary',num2str(bl-1),'.txt'],'w');
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
fclose(fid0);



