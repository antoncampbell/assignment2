function box_potential_2b
% clc;
% clear all;
% close all;

%Program Parameters
for multiplier=1:1:6
    uu=multiplier/4;
    
clear cMap
clear B
clear G


nx=200*uu;
ny=100*uu;
V0=1;

%Variables
sigma_box=1e-2;
sigma_main=1;
cMap=ones(ny,nx);
B=zeros(1,nx*ny);
G=sparse(nx*ny,nx*ny);


%Set conductivity map
cMap(:,:)=sigma_main;
cMap(1:40*uu,(80*uu+1):120*uu)=sigma_box;
cMap((60*uu+1):100*uu,(80*uu+1):120*uu)=sigma_box;

for i=1:nx
    for j=1:ny
        N=j+(i-1)*ny;
        
        if(i==1)%left
            
            %G(N,:)=0;
            G(N,N)=1;
            B(N)=V0;
            
        elseif(i==nx)%right
            
            %G(N,:)=0;
            G(N,N)=1;
            B(N)=0;
            
        elseif(j==1)%bottom
            Nyu=j+1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;
            
            Ryu=(cMap(j,i)+cMap(j+1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryu+Rxu+Rxd);
            G(N,Nyu)=Ryu;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
            
        elseif(j==ny)%top
            Nyd=j-1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;
            
            Ryd=(cMap(j,i)+cMap(j-1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryd+Rxu+Rxd);
            G(N,Nyd)=Ryd;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*ny;
            Nyd=j-1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;

            Ryu=(cMap(j,i)+cMap(j+1,i))/2.0;
            Ryd=(cMap(j,i)+cMap(j-1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryu+Ryd+Rxu+Rxd);
            G(N,Nyu)=Ryu;
            G(N,Nyd)=Ryd;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
        end
            
            
            
        
    end
end



V=G\B';

for i=1:nx
    for j=1:ny


        N=j+(i-1)*ny;

        Vmat(j,i)=V(N);

    end
end

[Ex,Ey]=gradient(-Vmat);
Emag=sqrt(Ex.^2+Ey.^2);
CurrentLeft(multiplier)=sum(cMap(:,1).*Ex(:,1));
CurrentRight(multiplier)=sum(cMap(:,nx).*Ex(:,nx));


% startlabel=3;
% figure(startlabel+1);
% surf(cMap,'EdgeColor','flat');
% title('Conductivity Map');
% ylabel('y position');
% xlabel('x position');
% cb4=colorbar;
% cb4.Label.String = '\sigma (S/m)';
% view(0,90);




pause(1)
%uu

end

figure(8);
hold on;
plot(CurrentLeft,'LineStyle','-','Marker','+','Color','blue');
plot(CurrentRight,'LineStyle','-.','Marker','o','Color','red');
legend('Left Contact','Right Contact');
hold on;
title('Figure 8: Current vs. Mesh Size');
ylabel('Relative Current');
xlabel('Relative Mesh Size Multiplier');

end

