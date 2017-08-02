function [Wp,sol,sys] = InitWFSim(Wp,options)
    WFSim_addpaths; % Add paths

    % Create empty structs
    sys = struct;
    sol = struct;

    % Initialize time
    sol.k = 0;
    
    % Create meshing and import control settings
    [Wp] = meshing(Wp,options.plotMesh,options.plotMesh); 

    % Initial flow fields
    [sol.u,sol.uu] = deal( Wp.site.u_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.v,sol.vv] = deal( Wp.site.v_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.p,sol.pp] = deal( Wp.site.p_init * ones(Wp.mesh.Nx,Wp.mesh.Ny) ); 

    if options.Linearversion
        sol.ul = sol.u;
        sol.vl = sol.v;
        sol.pl = sol.p;
        [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
    end;

    % Initialize parameters are empty matrices
    [sol.power,sol.ct,sol.Ueffect,sol.a] = deal(zeros(Wp.turbine.N,Wp.sim.NN)); 

    % Compute boundary conditions and matrices B1, B2
    [sys.B1,sys.B2,sys.bc]  = Compute_B1_B2_bc(Wp);
    sys.pRCM                = []; % Load empty vector
    
    % Compute projection matrices Qsp and Bsp
    if options.Projection
        [sys.Qsp, sys.Bsp]  = Solution_space(B1,B2,bc); % Projection matrices
        Wp.Nalpha           = size(sys.Qsp,2);
    end
end

