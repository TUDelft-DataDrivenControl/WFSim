function derivatives=derivatives_calc(A,Al,Ax,Ay,Wp,B1,b,bl,dSm,ccx,ccy,dbcdx) 

%% Remove columns and row and construct derivatives for optimization
    A(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    b(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    A(:,end) = [];A(end,:) = [];b(end)=[];
    A(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    
    Al(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    bl(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    Al(:,end) = [];Al(end,:) = [];bl(end)=[];
    Al(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    
    dAdx     = Al;
    dSmdx = sparse((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2));
    dSmdx(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)) =[dSm.xdu dSm.xdv; dSm.ydu dSm.ydv];
    
%     dSmdx = sparse((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2));
%     
%     dSmdx(1:(Wp.Nx-3)*(Wp.Ny-2),1:2*(Wp.Nx-3)*(Wp.Ny-2)) = [dSm.xxdu dSm.xxdv];
%     dSmdx((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),...
%         1:2*(Wp.Nx-2)*(Wp.Ny-3)) = [dSm.yydu dSm.yydv];
    dSmdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    dSmdx(:,end) = [];
    dSmdx(end,:) = [];
    dSmdx(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    
    
    dSm.dJdx(end,:)=[];
    dSm.dJdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:)=[];
    
    dSmdbeta = sparse((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2),Wp.N);
    dSmdbeta(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),:) = dSm.dbeta;
    dSmdbeta(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    dSmdbeta(end,:) = [];
    
    dSmdphi = sparse((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2),Wp.N);
    dSmdphi(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),:) = dSm.dphi;
    dSmdphi(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    dSmdphi(end,:) = [];
    
        dSmdPower_in = sparse((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2),Wp.N);
    dSmdPower_in(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),:) = dSm.dbeta.*dSm.dPower_in;
    dSmdPower_in(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    dSmdPower_in(end,:) = [];
    
    %dSmd  = blkdiag(dSm.uudx,dSm.vvdx,sparse((Wp.Ny-2)*(Wp.Nx-2)-2,(Wp.Ny-2)*(Wp.Nx-2)-2));
    Q     = blkdiag(diag(ccx),diag(ccy),sparse((Wp.Ny-2)*(Wp.Nx-2)-2,(Wp.Ny-2)*(Wp.Nx-2)-2));
    
    dbcdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    dbcdx(:,end) = [];dbcdx(end,:) = [];
    dbcdx(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    %% Store derivatives
    derivatives.dAdx        = dAdx;
    derivatives.dSmdbeta    = dSmdbeta;
    derivatives.dSmdphi     = dSmdphi;
       derivatives.dSmdPower_in     = dSmdPower_in;
    derivatives.dSmdx       = dSmdx;
    derivatives.Q           = Q;
    derivatives.dBc         = dbcdx;
    derivatives.dJdx        = dSm.dJdx;
    derivatives.dJdbeta     = dSm.dJdbeta;
    derivatives.dJdPhi     = dSm.dJdPhi;
    derivatives.dJdPower_in     = dSm.dJdPower_in;
    derivatives.A           = A;