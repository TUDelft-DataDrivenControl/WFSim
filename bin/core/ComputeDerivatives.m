function derivatives = ComputeDerivatives(A,Al,dSm,StrucDynamical,StrucBCs,Wp)

Ny    = Wp.mesh.Ny;
Nu    = Wp.Nu;
Nv    = Wp.Nv;
Np    = Wp.Np;

N     = Wp.turbine.N;

ccx   = StrucDynamical.ccx;
ccy   = StrucDynamical.ccy;
dbcdx = StrucBCs.dbcdx; 

dAdx  = Al;

dSmdx = sparse(Nu+Nv+Np+2,Nu+Nv+Np+2);
dSmdx(1:Nu+Nv,1:Nu+Nv) =[dSm.xdu dSm.xdv; dSm.ydu dSm.ydv];

%% New uncomment
dSmdx(Nu+Nv+Np+2-(Ny-2)+1,:) = [];
dSmdx(:,end) = [];
dSmdx(end,:) = [];
dSmdx(:,Nu+Nv+Np+2-(Ny-2)+1) = [];

dSm.dJdx(end,:)=[];
dSm.dJdx(Nu+Nv+Np+2-(Ny-2)+1,:)=[];
%%

dSmdbeta = sparse(Nu+Nv+Np+2,N);
dSmdbeta(1:Nu+Nv,:) = dSm.dbeta;
%% New uncomment
dSmdbeta(Nu+Nv+Np+2-(Ny-2)+1,:) = [];
dSmdbeta(end,:) = [];
%%

dSmdphi = sparse(Nu+Nv+Np+2,N);
dSmdphi(1:Nu+Nv,:) = dSm.dphi;
%% New uncomment
dSmdphi(Nu+Nv+Np+2-(Ny-2)+1,:) = [];
dSmdphi(end,:) = [];
%%

Q = blkdiag(diag(ccx),diag(ccy),sparse(Np,Np));

%% New uncomment
dbcdx(Nu+Nv+Np+2-(Ny-2)+1,:) = [];
dbcdx(:,end) = [];dbcdx(end,:) = [];
dbcdx(:,Nu+Nv+Np+2-(Ny-2)+1) = [];

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
dSmdPower_in            = sparse(Nu+Nv+Np+2,N);
dSmdPower_in(1:Nu+Nv,:) = dSm.dbeta.*dSm.dPower_in;

dSmdPower_in(Nu+Nv+Np+2-(Ny-2)+1,:) = [];
dSmdPower_in(end,:) = [];

derivatives.dSmdPower_in = dSmdPower_in;
derivatives.dJdPower_in  = dSm.dJdPower_in;