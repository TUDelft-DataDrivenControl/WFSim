function [Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh)

sys    = struct;
sol    = struct;

Projection    = options.Projection;
Linearversion = options.Linearversion;

% Create meshing and import control settings
[Wp,input]   = meshing(Wp,plotMesh,1); 

% Initial flow fields
[sol.u,sol.uu] = deal(Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
[sol.v,sol.vv] = deal(Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
[sol.p,sol.pp] = deal(Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny)); 


if Linearversion
    sol.ul = sol.u;
    sol.vl = sol.v;
    sol.pl = sol.p;
    [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
end;

% Initialize parameters are empty matrices
[Power,CT,Ueffect,a] = deal(zeros(Wp.turbine.N,Wp.sim.NN)); 

% Compute boundary conditions and matrices B1, B2
[B1,B2,bc]           = Compute_B1_B2_bc(Wp);
B2                   = 2*B2;

% Compute projection matrices Qsp and Bsp
if Projection
    [sys.Qsp, sys.Bsp]  = Solution_space(B1,B2,bc); % Projection matrices
    Wp.Nalpha           = size(sys.Qsp,2);
end;

end

