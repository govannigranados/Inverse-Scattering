clear
%%%    Copyright (C) 2022  Isaac Harris + Govanni Granados

tic
%% These are the keys to the car
centers = [0 0.5; 0 -0.75; -0.25 -0.25; 0.25 0.25]; %each row contains the centers of the defective region
rho = [0.95 1 0.9 1.1]; %used in functions

% each rho is the forcing term for the corresping center region 
% length(rho) must equal #of centers

%% Set the temperature

xs=linspace(-1,1,100); %meshgrid
Nang=64; %number of angles around unit circle
ka=25; %wave number
epsilon = 0.01; %radius of subregions

ang=linspace(-pi,pi,Nang+1); 
ang=ang(1:Nang);
yhat = [cos(ang); sin(ang)]; % yhat points on the boundary

% function takes in ingredients from above
[u, du] = mdata2(centers,ka,epsilon,yhat,rho); 

%%%%%%%%%%%%%%%% scatfield with error %%%%%%%%%%%%%%%
Eu=(-1 + 2.*rand(Nang,1)); Eu=(Eu+1i*Eu)/(sqrt(2)*norm(Eu)); 
Edu=(-1 + 2.*rand(Nang,1)); Edu=(Edu+1i*Edu)/(sqrt(2)*norm(Edu)); 
delta=0.25; U=u.*(1+delta*Eu); dU = du.*(1+delta*Edu);

%%%%%%%%%% 
pwave=@(z1,z2)  exp(1i*ka*(z1.*yhat(1,:)+z2.*yhat(2,:)));
dpwave=@(z1,z2) 1i*ka*(z1.*yhat(1,:)+z2.*yhat(2,:)).*exp(1i*ka.*(z1.*yhat(1,:)+z2.*yhat(2,:)));

Rh=@(z1,z2) dot(U,dpwave(z1,z2))-dot(dU,pwave(z1,z2));

Rhvec=zeros(Nang,1);

for ir=1:Nang
    Rhvec(ir)=Rh(yhat(1,ir),yhat(2,ir));
end 

Wh=zeros(max(size(xs)),max(size(xs)));


for i=1:max(size(xs))
    for j=1:max(size(xs)) 
        Wh(i,j)=abs(dot(Rhvec,(pwave(xs(i),xs(j))))); 
    end 
end 

% Normalization of imaging functional
Wh=Wh./(max(Wh(:)));
p=4;

%Outer Boundary given by yhat(1,:), yhat(2,:)

%%%%%% Graphic visualization 
figure;
subplot(1,2,1)
contour(xs,xs,(Wh').^p,'fill','on')
colormap hot
colorbar
hold on 
%plot(xc,yc,'--w','LineWidth',3) 
plot(yhat(1,:),yhat(2,:),'--w','LineWidth',3 )
axis('square')
hold off 

subplot(1,2,2)
surf(xs,xs,(Wh').^p)
%colorbar
shading interp
hold on 
plot(yhat(1,:),yhat(2,:),'--w','LineWidth',3) 
axis('square')
hold off 
toc