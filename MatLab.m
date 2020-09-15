clc
clear all
parpool (4)
rho=1;
Re=10;
Nx=128-2;
Ny=64-2;
dt=0.01;
tend=0.5;
L=2;
h=1;
dx=L/Nx;
dy=h/Ny;
init(1:Ny+2,1:Nx+2)=0;
uinit(1:Ny+2,1:Nx+2)=0;
B=dx/dy;
flag1=1;
flag=1;
x=dx/2:dx:L-dx/2;
y=dy/2:dy:h-dy/2;
count1=0;
count2=0;

uinit(1:Ny+2,1)=1;  % inlet
uinit(1,1:Nx+2)=-uinit(2,1:Nx+2);% lower wall
uinit(Ny+2,1:Nx+2)=-uinit(Ny+1,1:Nx+2);  % upper wall
uinit(1:Ny+2,Nx+1)=uinit(1:Ny+2,Nx); %  outlet
       
        
     spmd
     parts = [(Nx+2)/numlabs (Nx+2)/numlabs (Nx+2)/numlabs (Nx+2)/numlabs];
     numLocalCols = parts(labindex);
     leftColInd = sum(parts(1:labindex - 1)) + 1;
     rightColInd = leftColInd + numLocalCols - 1;
     u = uinit(:, leftColInd:rightColInd);
     v = init(:, leftColInd:rightColInd);
     Pp = init(:, leftColInd:rightColInd);
     P = init(:, leftColInd:rightColInd);
     if (labindex > 1),       u = [zeros(Ny+2, 1) u];  v = [zeros(Ny+2, 1) v]; ...
             Pp = [zeros(Ny+2, 1) Pp];  P = [zeros(Ny+2, 1) P]; end
     if (labindex < numlabs), u = [u zeros(Ny+2, 1)];  v = [v zeros(Ny+2, 1)]; ...
             Pp = [Pp zeros(Ny+2, 1)];  P = [P zeros(Ny+2, 1)]; end
     if (labindex == 1) || (labindex == numlabs)
    numLocalCols = numLocalCols - 1;
     end
       
     nooflabs=numlabs;
     up=u; vp=v; C=u;  % declaring array
     rightNeighbor = mod(labindex, numlabs) + 1;
     leftNeighbor  = mod(labindex - 2, numlabs) + 1;
    north    = 1:Ny;
    south    = 3:Ny + 2;
    currRow  = 2:Ny + 1;
    currCol  = 2:numLocalCols + 1;
    east     = 3:numLocalCols + 2;
    west     = 1:numLocalCols;
     end

tic
      
for t=0:dt:tend
    
     spmd
    t
    flag1 = 1;
    flag=1;
    %%
    while flag1 == 1  
        
        count2=count2+1;
    
 up(currRow,currCol)=u(currRow,currCol)-dt/dx/4*((u(currRow,currCol)+u(currRow,east)).^2-(u(currRow,currCol)+u(currRow,west)).^2)...
         -dt/dy/4*((u(currRow,currCol)+u(south,currCol)).*(v(currRow,currCol)+v(currRow,east))-(u(currRow,currCol)+u(north,currCol))...
         .*(v(north,currCol)+v(north,east)))-(P(currRow,east)-P(currRow,currCol))*dt/dx+dt/Re*((u(currRow,east)-2*u(currRow,currCol)+...
         u(currRow,west))/dx^2+(u(currRow+1,currCol)-2*u(currRow,currCol)+u(north,currCol))/dy^2);
 
     rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,up(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,up(:,end-1));
    if (labindex > 1),       up(:, 1) = leftBoundary;    end
    if (labindex < numlabs), up(:, end) = rightBoundary; end
    if labindex==1 , up(1:Ny+2,1)=1;  end
       up(1,:)=-up(2,:);
       up(Ny+2,:)=-up(Ny+1,:);
    if labindex==numlabs ,up(1:Ny+2,end)=up(1:Ny+2,end-1);   end
  
    
vp(currRow,currCol)=v(currRow,currCol)-dt/dy/4*((v(currRow,currCol)+v(south,currCol)).^2-(v(currRow,currCol)+v(north,currCol)).^2)...
        -dt/dx/4*((u(currRow,currCol)+u(south,currCol)).*(v(currRow,currCol)+v(currRow,east))-(u(currRow,west)+u(south,west))...
        .*(v(currRow,currCol)+v(currRow,west)))-(P(south,currCol)-P(currRow,currCol))*dt/dy+dt/Re*((v(currRow,east)-2*v(currRow,currCol)...
        +v(currRow,west))/dx^2+(v(south,currCol)-2*v(currRow,currCol)+v(north,currCol))/dy^2);
    
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,vp(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,vp(:,end-1));
    if (labindex > 1),       vp(:, 1) = leftBoundary;    end
    if (labindex < numlabs), vp(:, end) = rightBoundary; end

    if labindex==1 , vp(1:Ny+2,1)=-vp(1:Ny+2,2);  end
    if labindex==numlabs ,vp(1:Ny+2,end)=vp(1:Ny+2,end-1);   end   % I changed the indexs from the originial one
       vp(1,:)=0;% lower wall
       vp(Ny+1,:)=0;  % upper wall


    %LSOR
   C(currRow,currCol)=((up(currRow,currCol)-up(currRow,west))*dx/dt+(vp(currRow,currCol)-vp(north,currCol))*B^2*dy/dt);
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,C(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,C(:,end-1));
    if (labindex > 1),       C(:, 1) = leftBoundary;    end
    if (labindex < numlabs), C(:, end) = rightBoundary; end

  while flag == 1  
          count1=count1+1;
    Pp_temp=Pp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,Pp(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,Pp(:,end-1));
    if (labindex > 1),       Pp(:, 1) = leftBoundary;    end
    if (labindex < numlabs), Pp(:, end) = rightBoundary; end
    Pp(currRow,currCol)=0.25*(Pp(currRow,east)+Pp(currRow,west)+B^2*(Pp(south,currCol)+Pp(north,currCol))-C(currRow,currCol));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  consider change Pp to Pp_temp for velocity
    R_temp=abs((Pp-Pp_temp)./(Pp_temp+1e-9));
    Rw=max(max(R_temp));
    R=gplus(Rw);
    if R>=0.001 , flag=1;   else  flag=2;  end
  end
   flag=1;

    P_temp=P;
    P(currRow,currCol)=P(currRow,currCol)+Pp(currRow,currCol);
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,P(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,P(:,end-1));
    if (labindex > 1),       P(:, 1) = leftBoundary;    end
    if (labindex < numlabs), P(:, end) = rightBoundary; end
    
    

    R1_temp=abs((P-P_temp)./(P_temp+1e-9));
    R1w=max(max(R1_temp));
    R1=gplus(R1w);
    if R1>=0.00001, flag1=1;  else  flag1=2; end 
    
     end
   
v(currRow,currCol)=vp(currRow,currCol)+(Pp(currRow,currCol)-Pp(south,currCol))*dt/dy;
u(currRow,currCol)=up(currRow,currCol)+(Pp(currRow,currCol)-Pp(currRow,east))*dt/dx;
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,v(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,v(:,end-1));
    if (labindex > 1),       v(:, 1) = leftBoundary;    end
    if (labindex < numlabs), v(:, end) = rightBoundary; end
    rightBoundary = labSendReceive(leftNeighbor,rightNeighbor,u(:,2));
    leftBoundary = labSendReceive(rightNeighbor,leftNeighbor,u(:,end-1));
    if (labindex > 1),       u(:, 1) = leftBoundary;    end
    if (labindex < numlabs), u(:, end) = rightBoundary; end

    end   
     udraw=[up{:}];
     udraw= (udraw(2:end-1,2:end-1)+udraw(2:end-1,3:end))/2;
     surf(udraw)
     drawnow
end 
toc
delete(gcp)
