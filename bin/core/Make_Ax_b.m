function [sol, sys] = Make_Ax_b(Wp,sys,sol,options)
Nx    = Wp.mesh.Nx;
Ny    = Wp.mesh.Ny;

% Decide whether to start from uniform flow field or steady state
if sol.k == 1 && options.startUniform == 0
    dt = Inf;
else
    dt = Wp.sim.h;
end;

[StrucDiscretization]                = SpatialDiscr_Hybrid(Wp,sol,options.Linearversion); % Spatial discretization
[StrucDiscretization,StrucDynamical] = Dynamical(Wp,StrucDiscretization,sol,dt,options.Linearversion);  % Dynamical term
[StrucActuator,sol]                  = Actuator(Wp,sol,options); % Actuator/forcing function
[StrucDiscretization,StrucBCs]       = BoundaryConditions(Wp,StrucDiscretization,sol,options.Linearversion); % Zero gradient boundary conditions momentum equations

% Setup A matrix
Ay    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ay,2,3,1);
Ax    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ax,3,2,1);
sys.A = [blkdiag(Ax,Ay) [sys.B1;sys.B2]; ...
        [sys.B1;2*sys.B2]' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))]; 
sys.M = blkdiag(spdiags(StrucDynamical.ccx,0,Wp.Nu,Wp.Nu),...
                spdiags(StrucDynamical.ccy,0,Wp.Nv,Wp.Nv),...
                spdiags(zeros(Wp.Np,1),0,Wp.Np,Wp.Np));

if  options.Projection
    sys.Ct = blkdiag(spdiags(StrucDynamical.ccx,0,Wp.Nu,Wp.Nu),spdiags(StrucDynamical.ccy,0,Wp.Nv,Wp.Nv));
    sys.Et = sys.Qsp'*blkdiag(Ax,Ay)*sys.Qsp;
    sys.At = sys.Qsp'*sys.Ct*sys.Qsp;
    sys.St = sys.Qsp'*[StrucBCs.bx;StrucBCs.by] - sys.Qsp'*blkdiag(Ax,Ay)*sys.Bsp + sys.Qsp'*sys.Ct*sys.Bsp;
    sys.Bt = sys.Qsp'*[StrucActuator.Sm.xx;StrucActuator.Sm.yy];
    
    if k==1
        sys.pRCM = symrcm(sys.Et); % Calculate RCM
    end
    
    if options.Linearversion
        Ayl         = MakingSparseMatrixl( Nx,Ny,StrucDiscretization.day,2,3,1);
        Axl         = MakingSparseMatrixl( Nx,Ny,StrucDiscretization.dax,3,2,1);
        [Axlo,Aylo] = MakingSparseMatrixlo(Nx,Ny,StrucDiscretization.dax,StrucDiscretization.day);
        sys.Al      = [Axl Axlo sys.B1;Aylo Ayl sys.B2;sys.B1' 2*sys.B2' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))]-sys.A;
        
        sys.Atl = StrucDynamical.dcdx+StrucBCs.dbcdx+StrucActuator.dSm.dx-sys.Al;
        sys.Atl = sys.Qsp'*sys.Atl(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3))*sys.Qsp;
        sys.Etl = sys.Et;
        sys.Btl = sys.Qsp'*[StrucActuator.Sm.dxx;StrucActuator.Sm.dyy];
        % Model = dss(Atl,Btl,Qsp,0,Etl,1);
    end
    
else
    sys.b    = [StrucBCs.bx+StrucDynamical.cx+vec(StrucActuator.Sm.x');
        StrucBCs.by+StrucDynamical.cy+vec(StrucActuator.Sm.y');
        sys.bc];
    sys.m    = [StrucBCs.bx+vec(StrucActuator.Sm.x');
        StrucBCs.by+vec(StrucActuator.Sm.y');
        sys.bc];
    
    if options.Linearversion
        % Linear version
        Ayl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.day,2,3,1);
        Axl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.dax,3,2,1);
        [Axlo,Aylo] = MakingSparseMatrixlo(Nx,Ny,StrucDiscretization.dax,StrucDiscretization.day);
        Al          = [Axl Axlo sys.B1;Aylo Ayl sys.B2;sys.B1' 2*sys.B2' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))]-sys.A;
        
        sys.Al = (StrucDynamical.dcdx+StrucBCs.dbcdx+StrucActuator.dSm.dx-Al);
        sys.Bl = [StrucActuator.Sm.dxx;StrucActuator.Sm.dyy;zeros(length(sys.bc),size(input.beta,1)+size(input.phi,1))];
        sys.bl = [vec(StrucActuator.Sm.dx');
            vec(StrucActuator.Sm.dy');
            0.*sys.bc];
        
        sys.Al(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
        sys.bl(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
        sys.Bl(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
        sys.Al(:,size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1) = [];
        sys.Al(:,end) = [];sys.Al(end,:) = [];sys.bl(end) = [];sys.Bl(end,:) = [];
    end
    
    sys.A(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
    sys.b(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
    sys.m(size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1,:) = [];
    sys.A(:,size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1) = [];
    sys.A(:,end) = [];sys.A(end,:) = [];sys.b(end)=[];sys.m(end)=[];
    
    if sol.k==1
        sys.pRCM = symrcm(sys.A);
    end;
    
end

if (options.Linearversion>0 && options.Derivatives>0)
    sys.derivatives = ComputeDerivatives(sys.A,sys.Al,StrucActuator.dSm,StrucDynamical,StrucBCs,Wp) ;
end