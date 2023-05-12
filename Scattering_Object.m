classdef Scattering_Object
    %SCATTERING_OBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    %%%    Copyright (C) 2022  Isaac Harris + Govanni Granados
    centers
    rho
    xs
    Nang
    ka
    epsilon
    ang
    yhat
    u
    du
    Eu
    Edu
    delta
    pwave
    pdwave
    Rh
    Rhvec
    Wh
    p

    end
    
    methods
    function obj = Scattering_Object(~, ~)
            %SCATTERING_OBJECT Construct an instance of this class;
            obj.centers = [0 0.5; 0 -0.75; -0.25 -0.25; 0.25 0.25]; %each row contains the centers of the defective region
            rho = [0.95 1 0.9 1.1]; %used in functions
            %% Set the temperature
        
            obj.xs=linspace(-1,1,100); %meshgrid
            Nang=64; %number of angles around unit circle
            ka=25; %wave number
            epsilon = 0.01; %radius of subregions
            
            ang=linspace(-pi,pi,Nang+1); 
            ang=ang(1:Nang);
            obj.yhat = [cos(ang); sin(ang)]; % yhat points on the boundary
            
            % function takes in ingredients from above
            [u, du] = mdata2(obj.centers,ka,epsilon,obj.yhat,rho); 
            
            %%%%%%%%%%%%%%%% scatfield with error %%%%%%%%%%%%%%%
            Eu=(-1 + 2.*rand(Nang,1)); Eu=(Eu+1i*Eu)/(sqrt(2)*norm(Eu)); 
            Edu=(-1 + 2.*rand(Nang,1)); Edu=(Edu+1i*Edu)/(sqrt(2)*norm(Edu)); 
            delta=0.25; U=u.*(1+delta*Eu); dU = du.*(1+delta*Edu);
            
            %%%%%%%%%% 
            pwave=@(z1,z2)  exp(1i*ka*(z1.*obj.yhat(1,:)+z2.*obj.yhat(2,:)));
            dpwave=@(z1,z2) 1i*ka*(z1.*obj.yhat(1,:)+z2.*obj.yhat(2,:)).*exp(1i*ka.*(z1.*obj.yhat(1,:)+z2.*obj.yhat(2,:)));
            
            Rh=@(z1,z2) dot(U,dpwave(z1,z2))-dot(dU,pwave(z1,z2));
            
            Rhvec=zeros(Nang,1);
            
            for ir=1:Nang
                Rhvec(ir)=Rh(obj.yhat(1,ir),obj.yhat(2,ir));
            end 
            
            Wh=zeros(max(size(obj.xs)),max(size(obj.xs)));
            
            
            for i=1:max(size(obj.xs))
                for j=1:max(size(obj.xs)) 
                    obj.Wh(i,j)=abs(dot(Rhvec,(pwave(obj.xs(i),obj.xs(j))))); 
                end 
            end 
            
            % Normalization of imaging functional
            obj.Wh=obj.Wh./(max(obj.Wh(:)));
            obj.p=4;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        function Visualize(obj)
        %%%%%% Graphic visualization 
        figure;
        subplot(1,2,1)
        contour(obj.xs,obj.xs,(obj.Wh').^obj.p,'fill','on')
        colormap hot
        colorbar
        shading interp
        hold on 
        %plot(xc,yc,'--w','LineWidth',3) 
        plot(obj.yhat(1,:),obj.yhat(2,:),'--w','LineWidth',3 )
        axis('square')
        hold off 

        subplot(1,2,2)
        surf(obj.xs,obj.xs,(obj.Wh').^obj.p)
        %colorbar
        shading interp
        hold on 
        plot(obj.yhat(1,:),obj.yhat(2,:),'--w','LineWidth',3) 
        axis('square')
        hold off 
        end
    end
end

