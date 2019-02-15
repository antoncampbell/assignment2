function box_potential_2c
% clc;
% clear all;
% close all;

%Program Parameters
nx=200;
ny=100;
V0=1;

%Variables
sigma_box=1e-2;
sigma_main=1;
cMap=ones(ny,nx);
B=zeros(1,nx*ny);
G=sparse(nx*ny,nx*ny);

MaxLbDiv2=25;
MaxWb=49;
multiplier=2;
VaryWb_CurrentLeft=zeros(1,MaxWb);
VaryWb_CurrentRight=zeros(1,MaxWb);
VaryLb_CurrentLeft=zeros(1,MaxLbDiv2);
VaryLb_CurrentRight=zeros(1,MaxLbDiv2);

LbD2=20;
for Wb=1:MaxWb
       
        
%Set conductivity map
cMap(:,:)=sigma_main;
cMap(1:Wb,(nx/2-LbD2+1):(nx/2+LbD2))=sigma_box;
cMap((100-Wb+1):100,(nx/2-LbD2+1):(nx/2+LbD2))=sigma_box;

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



VaryWb_CurrentLeft(Wb)=sum(cMap(:,1).*Ex(:,1));
VaryWb_CurrentRight(Wb)=sum(cMap(:,nx).*Ex(:,nx));

end


figure(9)
hold on;
plot(VaryWb_CurrentLeft,'LineStyle','-','Color','blue');
plot(VaryWb_CurrentRight,'LineStyle','-.','Color','red');
legend('Left Contact','Right Contact');
hold off;
title('Figure 9: Current vs. Box Heights')
xlabel('Wb - Height of Each Box')
ylabel('Relative Current')




Wb=40;
for LbD2=1:1:MaxLbDiv2

        
        
%Set conductivity map
cMap(:,:)=sigma_main;
cMap(1:Wb,(nx/2-LbD2+1):(nx/2+LbD2))=sigma_box;
cMap((100-Wb+1):100,(nx/2-LbD2+1):(nx/2+LbD2))=sigma_box;

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


VaryLb_CurrentLeft(LbD2)=sum(cMap(:,1).*Ex(:,1));
VaryLb_CurrentRight(LbD2)=sum(cMap(:,nx).*Ex(:,nx));

end


figure(10)
hold on;
plot(2:2:2*MaxLbDiv2,VaryLb_CurrentLeft,'LineStyle','-','Color','blue');
plot(2:2:2*MaxLbDiv2,VaryLb_CurrentRight,'LineStyle','-.','Color','red');
legend('Left Contact','Right Contact');
hold off;
title('Figure 10: Current vs. Box Widths')
xlabel('Lb - Width of Each Box')
ylabel('Relative Current')

end





