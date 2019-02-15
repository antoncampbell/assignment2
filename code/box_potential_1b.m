function box_potential_1b

% clc;
% clear all;
% close all;

%Program Parameters
nx=50;
ny=60;
V0=1;

%Variables
cMap=ones(ny,nx);
B=zeros(1,nx*ny);
G=sparse(nx*ny,nx*ny);

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
            B(N)=V0;
            
        elseif(j==1)%bottom
            
            
%             Nyu=j+1+(i-1)*ny;
%             Nxu=j+(i-1+1)*ny;
%             Nxd=j+(i-1-1)*ny;
%             
%             Ryu=(cMap(j,i)+cMap(j+1,i))/2.0;
%             Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
%             Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
%             
%             %G(N,:)=0;
%             G(N,N)=-(Ryu+Rxu+Rxd);
%             G(N,Nyu)=Ryu;
%             G(N,Nxu)=Rxu;
%             G(N,Nxd)=Rxd;

            %G(N,:)=0;
            G(N,N)=1;
            B(N)=0;
            
        elseif(j==ny)%top
%             Nyd=j-1+(i-1)*ny;
%             Nxu=j+(i-1+1)*ny;
%             Nxd=j+(i-1-1)*ny;
%             
%             Ryd=(cMap(j,i)+cMap(j-1,i))/2.0;
%             Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
%             Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
%             
%             %G(N,:)=0;
%             G(N,N)=-(Ryd+Rxu+Rxd);
%             G(N,Nyd)=Ryd;
%             G(N,Nxu)=Rxu;
%             G(N,Nxd)=Rxd;

            %G(N,:)=0;
            G(N,N)=1;
            B(N)=0;
            
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

figure(2)
surf(Vmat);
title('Figure 2: Voltage Map for Finite Difference Method - Simple Case 2');
ylabel('y position');
xlabel('x position');
cb3=colorbar;
cb3.Label.String = 'V (V)';
view(0,90);

aa=ny-1;
bb=(nx-1)/2.0;
Ntot=116;
%nvec=1:1:100;
yvec=0:1:aa;
xvec=-bb:1:bb;
%[Nmesh,Xmesh,Ymesh]=meshgrid(nvec,xvec,yvec);
[Xmesh,Ymesh]=meshgrid(xvec,yvec);

Vanly=zeros(length(yvec),length(xvec));

for nn=1:2:Ntot
    
Vanly=Vanly+(4.*V0./pi).*(cosh(nn.*pi.*Xmesh./aa).*sin(nn.*pi.*Ymesh./aa)./(nn.*cosh(nn.*pi.*bb./aa)));

% figure(3)
% surf(Vanly)



end
figure(3)
surf(Vanly)
title('Figure 3: Voltage Map for Analytical Solution - Simple Case 2');
ylabel('y position');
xlabel('x position');
cb3=colorbar;
cb3.Label.String = 'V (V)';
view(0,90);

end
