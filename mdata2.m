%%% Code by Govanni Granados
% centers: contains the centers of the subregions in a matrix
% k: wave number
% epsilon: subregion radius
% yhat: [cos(ang) sin(ang)]

function [g_u, g_du] = mdata2(centers,k,epsilon,yhat,rho)

N = length(yhat); %= number of angles

us_integrand = @(p,x1,x2,y1,y2) -(1i/4).*p.*besselh(0,k.*sqrt((x1-y1).^2 + (x2 - y2).^2)); %changed - to + inside sqrt
nor_usint = @(p,x1,x2,y1,y2) (1i/4).*p.*k.*besselh(1,k.*sqrt((x1-y1).^2 + (x2 - y2).^2)).*((1-(x1.*y1 + x2.*y2))./(sqrt((x1-y1).^2 + (x2 - y2).^2)));

g_u = zeros(N,1);
g_du = zeros(N,1);

for n = 1:N
    s1 = 0;
    s2 = 0;
    for j = 1:size(centers,1)
        s1 = s1 + numquad2d(@(r,theta) us_integrand(rho(j),yhat(1,n),yhat(2,n),r.*cos(theta)+centers(j,1),r.*sin(theta)+centers(j,2)).*r,0,epsilon,-pi,pi);
        s2 = s2 + numquad2d(@(r,theta) nor_usint(rho(j),yhat(1,n),yhat(2,n),r.*cos(theta)+centers(j,1),r.*sin(theta)+centers(j,2)).*r,0,epsilon,-pi,pi);
    end
    g_u(n) = s1;
    g_du(n) = s2;
end