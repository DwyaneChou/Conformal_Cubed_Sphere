clc
clear

maxit = 100000;
Ermax = 1.E-15;
show  = 1;      % 1 for yes o for no  to disp solution while solving

dx        = 500000;
dlambda   = 1;
dtheta    = 1;

projection= 'Equiangular'; % Choose from 'Equidistant' or 'Equiangular'

R         = 1;
d2r       = pi/180;
r2d       = 180/pi;

dlambda   = dlambda*d2r;
dtheta    = dtheta *d2r;

a         = R/sqrt(3); % Length of the cubic edges

if strcmp(projection,'Equidistant')
    % Equidistant
    x         = -a:dx:a;
    y         = x;
elseif strcmp(projection,'Equiangular')
    % Equiangular
    lambda      = -pi/4:dlambda:pi/4;
    theta       = -pi/4 :dtheta :pi/4;
    x           = a*tan(lambda);
    y           = a*tan(theta );
end

nx        = size(x,2);
ny        = size(y,2);
a_matrix  = ones(nx,ny)*a;

x         = repmat(x,ny,1);
y         = repmat(y,nx,1)';

r         = sqrt(a_matrix.^2 + x.^2 + y.^2);

cart_coord(1,:,:) = R./r.*a_matrix;
cart_coord(2,:,:) = R./r.*x;
cart_coord(3,:,:) = R./r.*y;

% face 5
X = -squeeze(cart_coord(3,:,:));
Y = squeeze(cart_coord(2,:,:));
Z = squeeze(cart_coord(1,:,:));

% figure
% plot3(X,Y,Z,'.','Color','g');
% surf(X,Y,Z,'EdgeColor','k','FaceColor','r')
% hold on

newX = X;
newY = Y;
newZ = Z;
Er1  = zeros(1,maxit);
Er2  = zeros(1,maxit);
Er3  = zeros(1,maxit);
Err  = zeros(1,maxit);
pNorm = ones(size(X));
for t = 1:maxit
    i = 2:nx-1;
    j = 2:ny-1;
    
    newX(i,j) = ( X(i+1,j) + X(i-1,j) + X(i,j+1) + X(i,j-1) ) / 4;
    newY(i,j) = ( Y(i+1,j) + Y(i-1,j) + Y(i,j+1) + Y(i,j-1) ) / 4;
    newZ(i,j) = ( Z(i+1,j) + Z(i-1,j) + X(i,j+1) + Z(i,j-1) ) / 4;
    
    pNorm(i,j) = sqrt(newX(i,j).^2+newY(i,j).^2+newZ(i,j).^2);
    
    newX(i,j) = newX(i,j)./pNorm(i,j);
    newY(i,j) = newY(i,j)./pNorm(i,j);
    newZ(i,j) = newZ(i,j)./pNorm(i,j);
    
    Er1(1,t) = max(max(abs(newX-X)));
    Er2(1,t) = max(max(abs(newY-Y)));
    Er3(1,t) = max(max(abs(newZ-Z)));
    Err(1,t) = max([Er1(1,t),Er2(1,t),Er3(1,t)]);
    
    disp(['iter err=',num2str(Err(1,t))])
    
    X = newX;
    Y = newY;
    Z = newZ;
    
    if Err(t)<Ermax
        disp(['Convergence in ',num2str(Ermax),' by iter ',num2str(t)])
        break
    end
    
    if show==1
        if ceil(t/10)*10==t
            clf
            hold on
            axis equal
            title(['iter ',num2str(t),' err=',num2str(Err(t),'%4.2e')])
%             for m=1:nx
%                 plot3(X(m,:),Y(m,:),Z(m,:),'b');
%             end
%             for m=1:ny
%                 plot3(X(:,m),Y(:,m),Z(m,:),'Color',[0 0 0]);
%             end
            surf(X,Y,Z,'EdgeColor','k','FaceColor','r')
            pause(0.001)
        end
    end
end


