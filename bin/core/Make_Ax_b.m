function [sol, sys] = Make_Ax_b(Wp,sys,sol,options)
% Import variables
Nx    = Wp.mesh.Nx;
Ny    = Wp.mesh.Ny;

% Decide whether to start from uniform flow field or steady state
if sol.k == 1 && Wp.sim.startUniform == 0
    dt = Inf;
else
    dt = Wp.sim.h;
    dt = dt/2; % correction factor for difference with LES wake propagation
end;

% The following is really the core of WFSim: creating the system matrices
[StrucDiscretization]                = SpatialDiscr_Hybrid(Wp,sol,options.Linearversion); % Spatial discretization
[StrucDiscretization,StrucDynamical] = Dynamical(Wp,StrucDiscretization,sol,dt,options.Linearversion);  % Dynamical term
[StrucActuator,sol]                  = Actuator(Wp,sol,options); % Actuator/forcing function
[StrucDiscretization,StrucBCs]       = BoundaryConditions(Wp,StrucDiscretization,sol,options.Linearversion); % Zero gradient boundary conditions momentum equations

% Collect all terms and create the A matrix in 'A*x = b'
Ay    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ay,2,3,1);
Ax    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ax,3,2,1);


sys.A = [blkdiag(Ax,Ay) [sys.B1;sys.B2]; [sys.B1;2*sys.B2]' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))];       

% If necessary, project away the continuity equation in A*x = b            
if  options.Projection
    sys.Ct = blkdiag(spdiags(StrucDynamical.ccx,0,Wp.Nu,Wp.Nu),spdiags(StrucDynamical.ccy,0,Wp.Nv,Wp.Nv));
    sys.Et = sys.Qsp'*blkdiag(Ax,Ay)*sys.Qsp;
    sys.At = sys.Qsp'*sys.Ct*sys.Qsp;
    sys.St = sys.Qsp'*[StrucBCs.bx;StrucBCs.by] - sys.Qsp'*blkdiag(Ax,Ay)*sys.Bsp + sys.Qsp'*sys.Ct*sys.Bsp;
    sys.Bt = sys.Qsp'*[StrucActuator.Sm.xx;StrucActuator.Sm.yy];
    
    if sol.k==1
        sys.pRCM = symrcm(sys.Et); % Calculate RCM
    end
    
    if options.Linearversion
        Ayl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.day,2,3,1);
        Axl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.dax,3,2,1);
        [Axlo,Aylo] = MakingSparseMatrixlo(Nx,Ny,StrucDiscretization.dax,StrucDiscretization.day);
        sys.Al      = [Axl Axlo sys.B1;Aylo Ayl sys.B2;sys.B1' 2*sys.B2' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))]-sys.A;
        
        sys.Atl = StrucDynamical.dcdx+StrucBCs.dbcdx+StrucActuator.dSm.dx-sys.Al;
        sys.Atl = sys.Qsp'*sys.Atl(1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3),1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3))*sys.Qsp;
        sys.Etl = sys.Et;
        sys.Btl = sys.Qsp'*[StrucActuator.Sm.dxx;StrucActuator.Sm.dyy];
    end
    
else
    % Collect all terms and create the b vector in 'A*x = b'
    sys.b    = [StrucBCs.bx+StrucDynamical.cx+vec(StrucActuator.Sm.x');
                StrucBCs.by+StrucDynamical.cy+vec(StrucActuator.Sm.y');
                sys.bc];

    
    if options.Linearversion
        % Linear version
        Ayl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.day,2,3,1);
        Axl         = MakingSparseMatrixl(Nx,Ny,StrucDiscretization.dax,3,2,1);
        [Axlo,Aylo] = MakingSparseMatrixlo(Nx,Ny,StrucDiscretization.dax,StrucDiscretization.day);
        Al          = [Axl Axlo sys.B1;Aylo Ayl sys.B2;sys.B1' 2*sys.B2' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))]-sys.A;
        
        sys.Al = (StrucDynamical.dcdx+StrucBCs.dbcdx+StrucActuator.dSm.dx-Al);
        sys.Bl = [StrucActuator.Sm.dxx;StrucActuator.Sm.dyy;zeros(length(sys.bc),Wp.turbine.N*2)];
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
    sys.A(:,size(Ax,1)+size(Ay,1)+size(sys.B1',1)-(Ny-2)+1) = [];
    sys.A(:,end) = [];sys.A(end,:) = [];sys.b(end)=[];
    
    if sol.k==1
        sys.pRCM = symrcm(sys.A);
    end;
    
end
