function derivatives = ComputeDerivatives(A,Al,dSm,StrucDynamical,StrucBCs,Wp)

Nx    = Wp.mesh.Nx;
Ny    = Wp.mesh.Ny;
N     = Wp.turbine.N;

ccx   = StrucDynamical.ccx;
ccy   = StrucDynamical.ccy;
dbcdx = StrucBCs.dbcdx; 

dAdx  = Al;

dSmdx = sparse((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));
dSmdx(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)) =[dSm.xdu dSm.xdv; dSm.ydu dSm.ydv];

%dSmdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
%dSmdx(:,end) = [];
%dSmdx(end,:) = [];
%dSmdx(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1) = [];

%dSm.dJdx(end,:)=[];
%dSm.dJdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:)=[];

dSmdbeta = sparse((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2),N);
dSmdbeta(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),:) = dSm.dbeta;
%dSmdbeta(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
%dSmdbeta(end,:) = [];

dSmdphi = sparse((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2),N);
dSmdphi(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),:) = dSm.dphi;
%dSmdphi(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
%dSmdphi(end,:) = [];

Q = blkdiag(diag(ccx),diag(ccy),sparse((Ny-2)*(Nx-2)-2,(Ny-2)*(Nx-2)));

%dbcdx(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
%dbcdx(:,end) = [];dbcdx(end,:) = [];
%dbcdx(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1) = [];

%% Store derivatives
derivatives.ccx          = ccx;
derivatives.ccy          = ccy;
derivatives.dAdx         = dAdx;
derivatives.dSmdbeta     = dSmdbeta;
derivatives.dSmdphi      = dSmdphi;
derivatives.dSmdx        = dSmdx;
derivatives.Q            = Q;
derivatives.dBc          = dbcdx;
derivatives.dJdx         = dSm.dJdx;
derivatives.dJdbeta      = dSm.dJdbeta;
derivatives.dJdPhi       = dSm.dJdPhi;
derivatives.A            = A;


%% Before all is in
dSmdPower_in                                  = sparse((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2),N);
dSmdPower_in(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),:) = dSm.dbeta.*dSm.dPower_in;

%dSmdPower_in(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
%dSmdPower_in(end,:) = [];

derivatives.dSmdPower_in = dSmdPower_in;
derivatives.dJdPower_in  = dSm.dJdPower_in;