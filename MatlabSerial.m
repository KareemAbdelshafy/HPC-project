clc
clear all
%matlabpool (4)
rho=1;
Re=10;
Nx=20;
Ny=10;
dt=0.001
tend=0.5;
L=2;
h=1;
dx=L/Nx;
dy=h/Ny;
C(Ny+2,Nx+2)=0;
P(Ny+2,1:Nx+2)=4;
Po(1:Ny+2,1:Nx+2)=0;
Pa(1:Ny+2,1:Nx+2)=0;
Pb(1:Ny+2,1:Nx+2)=0;
u(1:Ny+2,1:Nx+2)=0;
ua(1:Ny+2,1:Nx+2)=0;
uo(1:Ny+2,1:Nx+2)=0;
v(1:Ny+2,1:Nx+2)=0;
va(1:Ny+2,1:Nx+2)=0;
vo(1:Ny+2,1:Nx+2)=0;
B=dx/dy;
w=1;
flag1=1;
flag=1;
x=dx/2:dx:L-dx/2;
y=dy/2:dy:h-dy/2;
count1=0;
count = 0;
count2=0;
R=0;
R1=0;
R2=0;
R3=0;
uc(Ny,Nx)=0;
vc(Ny,Nx)=0;
tic
        u(1:Ny+2,1)=1;  % inlet
        u(1,1:Nx+2)=-u(2,1:Nx+2);% lower wall
        u(Ny+2,1:Nx+2)=-u(Ny+1,1:Nx+2);  % upper wall
        u(1:Ny+2,Nx+1)=u(1:Ny+2,Nx); %  outlet
        uo=u;
        v(2:Ny+1,1)=-v(2:Ny+1,2);  % inlet
        v(1,2:Nx+1)=0;% lower wall
        v(Ny+2,2:Nx+1)=0;  % upper wall
        v(2:Ny+1,Nx+1)=v(2:Ny+1,Nx); %  outlet
        vo=v;
        
%         Pa(2,2:Ny+1)=0;
        Pa(2:Ny+1,1)=Pa(2:Ny+1,2); %inlet
        Pa(1,2:Nx+1)=Pa(2,2:Nx+1); %lower wall
        Pa(Ny+2,2:Nx+1)=Pa(Ny+1,2:Nx+1); %upper wall
        Pa(2:Ny+1,Nx+2)=2*Pa(2:Ny+1,Nx+1)-Pa(2:Ny+1,Nx); %outlet
figure        
for t=0:dt:tend
    t
    flag1 = 1;
    flag=1;
    count1=count1+1; 

    while flag1 == 1
     count2 = count2 +1;

%using Pinitial to calculate u and v
    %u(2:Nx,2:Ny)=uo(2:Nx,2:Ny);
    %v(2:Nx,2:Ny)=vo(2:Nx,2:Ny);
    
      for i=2:Nx+1
        for j=2:Ny+1
            uo(j,i)=u(j,i)-dt/dx/4*((u(j,i)+u(j,i+1))^2-(u(j,i)+u(j,i-1))^2)...
                -dt/dy/4*((u(j,i)+u(j+1,i))*(v(j,i)+v(j,i+1))-(u(j,i)+u(j-1,i))*(v(j-1,i)+v(j-1,i+1)))...
                -(P(j,i+1)-P(j,i))*dt/dx+dt/Re*((u(j,i+1)-2*u(j,i)+u(j,i-1))/dx^2+(u(j+1,i)-2*u(j,i)+u(j-1,i))/dy^2);
        end 
      end 
        uo(1:Ny+2,1)=1;
        uo(1,1:Nx+2)=-uo(2,1:Nx+2);
        uo(Ny+2,1:Nx+2)=-uo(Ny+1,1:Nx+2);
        uo(1:Ny+2,Nx+2)=uo(1:Ny+2,Nx+1);
                
      for i=2:Nx+1
        for j=2:Ny+1
            vo(j,i)=v(j,i)-dt/dy/4*((v(j,i)+v(j+1,i))^2-(v(j,i)+v(j-1,i))^2)...
                -dt/dx/4*((u(j,i)+u(j+1,i))*(v(j,i)+v(j,i+1))-(u(j,i-1)+u(j+1,i-1))*(v(j,i)+v(j,i-1)))...
                -(P(j+1,i)-P(j,i))*dt/dy+dt/Re*((v(j,i+1)-2*v(j,i)+v(j,i-1))/dx^2+(v(j+1,i)-2*v(j,i)+v(j-1,i))/dy^2);
        end 
      end
        vo(1:Ny+2,1)=-vo(1:Ny+2,2);
        vo(1,1:Nx+2)=0;% lower wall
        vo(Ny+1,1:Nx+2)=0;  % upper wall
        vo(1:Ny+2,Nx+2)=vo(1:Ny+2,Nx+1);
        
              C(2:Ny+1,2:Nx+1)=((uo(2:Ny+1,2:Nx+1)-uo(2:Ny+1,2-1:Nx+1-1))*dx/dt+(vo(2:Ny+1,2:Nx+1)-vo(1:Ny,2:Nx+1))*B^2*dy/dt);

    %LSOR
    while flag == 1 
    count= count+1;

        Pb=Pa;


                
        for i=2:Nx+1 
            for j= 2:Ny+1
                Pa(j,i)=0.25*(Pb(j,i+1)+Pb(j,i-1)+B^2*(Pb(j+1,i)+Pb(j-1,i))-C(j,i));
            end
             
        end
        
%         Pa(1,1:Nx+2)=0; 
%         Pa(Ny+2,1:Nx+2)=0;
%         Pa(2:Ny+1,1)=0;
%         Pa(2:Ny+1,Nx+2)=0;

     R_temp=abs((Pa-Pb)./Pb);
     R=max(max(R_temp));

    if R>=0.001 , flag=1;   else  flag=2;  end
    end 
    flag=1;

    
     Po=P;
     for i=2:Nx+1
         for j=2:Ny+1
            %vo(i,j)=vo(i,j)+va(i,j);
            %uo(i,j)=uo(i,j)+ua(i,j);
            P(j,i)=P(j,i)+Pa(j,i);
         end 
     end
     
     %R1=sqrt(sum(sum((P-Po).^2)))/sqrt(sum(sum((Po).^2)));
     %R2=sqrt(sum(sum(uo-u)).^2)/sqrt(sum(sum(u)).^2);
     %R3=sqrt(sum(sum(vo-v)).^2)/sqrt(sum(sum(v)).^2);
     
     R1_temp=abs((P-Po)./Po);
     R1=max(max(R1_temp));
     
    if R1>=0.00001
        flag1=1;
    else 
        flag1=2;
    end    
end 

       for i=2:Nx+1
           for j=2:Ny+1
               ua(j,i)=(Pa(j,i)-Pa(j,i+1))*dt/dx;
               va(j,i)=(Pa(j,i)-Pa(j+1,i))*dt/dy;
               v(j,i)=vo(j,i)+va(j,i);
               u(j,i)=uo(j,i)+ua(j,i);
           end 
       end
for i=2:Nx+1
    for j=2:Ny+1
        uc(j,i)=(uo(j,i)+uo(j,i+1))/2;
        vc(j,i)=(vo(j,i)+vo(j,i+1))/2;
    end 
end
surf(uc(2:Ny+1,2:Nx+1))
drawnow
end 
timeelapsed=toc
%matlabpool close

save serialS.mat

