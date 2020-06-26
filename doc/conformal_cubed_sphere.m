% clc
clear

iPatch = 2;

ccs_mesh = '..\mesh_generation\run\ccs_output.nc';
R2D = 180/pi;
D2R = pi/180;
x   = ncread(ccs_mesh,'x');
y   = ncread(ccs_mesh,'y');
z   = ncread(ccs_mesh,'z');
lon = ncread(ccs_mesh,'lon');
lat = ncread(ccs_mesh,'lat');

x1d   = reshape(x,[],1);
y1d   = reshape(y,[],1);
z1d   = reshape(z,[],1);
lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);

nx = size(x,1);
ny = size(x,2);
np = size(x,3);

radius = 1;
[xp,yp,zp] = sph2cart(lon1d*D2R,lat1d*D2R,radius);

% Plot white base
[xs,ys,zs] = sphere(800);
s = surf(xs*radius,ys*radius,zs*radius,'FaceColor','w');
set(s,'EdgeColor','None')
hold on

% scatter3(xp,yp,zp,'.')
scatter3(x1d,y1d,z1d,'.')

% figure
% plt = pcolor(squeeze(lat(:,:,iPatch)));
% set(plt,'EdgeColor','none')
% colormap(jet)

% patch(pointCoord(1,pointId_patch),pointCoord(2,pointId_patch),pointCoord(3,pointId_patch),areaCell(iCell));

% figure
% lon1d = reshape(lon(:,:,iPatch),[],1);
% lat1d = reshape(lat(:,:,iPatch),[],1);
% 
% [xp,yp,zp] = sph2cart(lon1d*D2R,lat1d*D2R,radius);
% scatter3(xp,yp,zp,'.')


% Calculate sphere mesh
xc = x(2:2:end,2:2:end,:);
yc = y(2:2:end,2:2:end,:);
zc = z(2:2:end,2:2:end,:);
nc = size(xc,1);
[xx,xy] = extend_mesh(xc);
[yx,yy] = extend_mesh(yc);
[zx,zy] = extend_mesh(zc);

dx  = 90/nc;
x1d = dx/2:90/nc:360-dx/2;
x1d = x1d * D2R;
dx  = dx*D2R;

for i = 1:nc
    pxx(i) = csape(x1d,xx(:,i),'periodic');
    pxy(i) = csape(x1d,xy(i,:),'periodic');
    pyx(i) = csape(x1d,yx(:,i),'periodic');
    pyy(i) = csape(x1d,yy(i,:),'periodic');
    pzx(i) = csape(x1d,zx(:,i),'periodic');
    pzy(i) = csape(x1d,zy(i,:),'periodic');

    dXdx(:,i) = pxx(i).coefs(1:nc,3);
    dYdx(:,i) = pyx(i).coefs(1:nc,3);
    dZdx(:,i) = pzx(i).coefs(1:nc,3);

    dXdy(i,:) = pxy(i).coefs(1:nc,3);
    dYdy(i,:) = pyy(i).coefs(1:nc,3);
    dZdy(i,:) = pzy(i).coefs(1:nc,3);
end

J1(:,:,1) = dXdx;
J1(:,:,2) = dYdx;
J1(:,:,3) = dZdx;

J2(:,:,1) = dXdy;
J2(:,:,2) = dYdy;
J2(:,:,3) = dZdy;

J = zeros(nc,nc,3,2);
G = zeros(nc,nc,2,2);
for j = 1:nc
    for i = 1:nc
        J(i,j,:,1) = J1(i,j,:);
        J(i,j,:,2) = J2(i,j,:);
        G(i,j,:,:) = squeeze(J(i,j,:,:))' * squeeze(J(i,j,:,:));
        sqrtG(i,j) = sqrt(det(squeeze(G(i,j,:,:))));
        area (i,j) = sqrtG(i,j) * dx^2;
    end
end

nuerical_sphere_area = 6 * sum(sum(area));
sphere_area_error    = 4*pi - nuerical_sphere_area;
disp(['Sphere area error = ',num2str(sphere_area_error,'%.15e')])

% Extend mesh
function [fx,fy] = extend_mesh(f)

nx = size(f,1);
ny = size(f,2);

fx = squeeze(f(:,:,1));
fy = squeeze(f(:,:,1));

% patch1
fx(1*nx+1:2*nx,:) = f(:,:,2);
fx(2*nx+1:3*nx,:) = f(:,:,3);
fx(3*nx+1:4*nx,:) = f(:,:,4);

fy(:,1*ny+1:2*ny) = f(:,:,5);
fy(:,2*ny+1:3*ny) = f(nx:-1:1,ny:-1:1,3);
fy(:,3*ny+1:4*ny) = f(:,:,6);

end