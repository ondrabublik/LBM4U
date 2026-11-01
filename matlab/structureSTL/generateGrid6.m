function generateGrid6

ni = 10;
nj = 10;
nk = 10;

name = 'forma2.stl';
[f,v] = stlread(name);

fv.faces = f;
fv.vertices = v;

figure
hold on;
patch(fv,'FaceColor','g','FaceAlpha',0.3,'linestyle','none')
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

n = 7;
dl = dx/n;
disp(['dl = ',num2str(dl)]);

n1 = n+1;
n2 = round(dy/dl) + 1;
n3 = round(dz/dl) + 1;

nc = (n1+1)*(n2+1)*(n3+1);
x = zeros(nc,1);
y = zeros(nc,1);
z = zeros(nc,1);

for i = 1:n1+1
    for j = 1:n2+1
        for k = 1:n3+1
            s = (n1+1)*(n2+1)*(k-1) + (n1+1)*(j-1) + i;
            x(s) = xmin - dl/4 + (i-1)*dl;
            y(s) = ymin - dl/4 + (j-1)*dl;
            z(s) = zmin - dl/4 + (k-1)*dl;
        end
    end
end

in = inpolyhedron(fv,[x,y,z]);
ex = [0 1 1 0 0 1 1 0];
ey = [0 0 1 1 0 0 1 1];
ez = [0 0 0 0 1 1 1 1];

% print blocks
hold on
nBlock = 0;
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            ind = [];
            b = 0;
            for p = 1:8
                indp = (n1+1)*(n2+1)*(k-1 + ez(p)) + (n1+1)*(j-1 + ey(p)) + i + ex(p);
                ind = [ind,indp];
                if(in(indp) == 1)
                    b = 1;
                end
            end
            if(b == 1)
                plot3(x([ind(1:4),ind(1)]),y([ind(1:4),ind(1)]),z([ind(1:4),ind(1)]),'k')
                plot3(x([ind(5:8),ind(5)]),y([ind(5:8),ind(5)]),z([ind(5:8),ind(5)]),'k')
                plot3(x([ind(1),ind(5)]),y([ind(1),ind(5)]),z([ind(1),ind(5)]),'k')
                plot3(x([ind(2),ind(6)]),y([ind(2),ind(6)]),z([ind(2),ind(6)]),'k')
                plot3(x([ind(3),ind(7)]),y([ind(3),ind(7)]),z([ind(3),ind(7)]),'k')
                plot3(x([ind(4),ind(8)]),y([ind(4),ind(8)]),z([ind(4),ind(8)]),'k')
                
                nBlock = nBlock+1;
            end
        end
    end
end
nBlock

mkdir('mesh');
mkdir('mesh/boundary');

disp('Buidling structure ...')
fid0 = fopen('mesh/blocks.txt','w');
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            b = 0;
            for p = 1:8
                ind = (n1+1)*(n2+1)*(k-1 + ez(p)) + (n1+1)*(j-1 + ey(p)) + i + ex(p);
                if(in(ind) == 1)
                    b = 1;
                end
            end
            if(b == 1)
                nib = ni;
                njb = nj;
                nkb = nk;
                
                nct = nib*njb*nkb;
                xt = zeros(nct,1);
                yt = zeros(nct,1);
                zt = zeros(nct,1);

                r = (n1+1)*(n2+1)*(k-1) + (n1+1)*(j-1) + i;
                for it = 1:nib
                    for jt = 1:njb
                        for kt = 1:nkb
                            s = nib*njb*(kt-1) + nib*(jt-1) + it;
                            xt(s) = x(r) + (it-1)/(nib-1)*dl + 1/(nib-1)*dl/2;
                            yt(s) = y(r) + (jt-1)/(njb-1)*dl + 1/(nib-1)*dl/2;
                            zt(s) = z(r) + (kt-1)/(nkb-1)*dl + 1/(nib-1)*dl/2;
                        end
                    end
                end
                int = inpolyhedron(fv,[xt,yt,zt]);
                plot3(xt(int),yt(int),zt(int),'.k')
                %plot3(xt(int==0),yt(int==0),zt(int==0),'.r')

                if(sum(int) > 0)
                    disp(['Block: ',num2str(n1*n2*(k-1) + n1*(j-1) + i-1)])
                    fprintf(fid0,'%d %d %d\n',i-1,j-1,k-1);
                    
                    bl = n1*n2*(k-1) + n1*(j-1) + i - 1;
                    fid = fopen(['mesh/boundary/blockBoundary',num2str(bl),'.txt'],'w');
                    for it = 1:nib
                        for jt = 1:njb
                            for kt = 1:nkb
                                b = 1;
                                ind = nib*njb*(kt-1) + nib*(jt-1) + it;
                                if(int(ind) == 1)
                                    b = 0;
                                end
                                fprintf(fid,'%i \n',b);
                            end
                        end
                    end
                    fclose(fid);
                end
            end
        end
    end
end
fclose(fid0);






