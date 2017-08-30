function [Wp,sol,sys] = InitWFSim(Wp,options)
%INITWFSIM  Initializes the WFSim model

    % Create empty structs
    sys = struct; % This struct will contain all the system matrices at time k
    sol = struct; % This struct will contain the solution (flowfields, power, ...) at time k
    
    % Import simulation scenario (meshing, atmospheric properties, turbine settings)
    [Wp] = meshing(Wp.name,options.plotMesh,options.plotMesh); 

    % Initialize time vector for sol at time k = 0
    sol = struct('k',0,'time',Wp.sim.time(1));
    
    % Initialize flow fields as uniform ('no turbines present yet')
    [sol.u,sol.uu] = deal( Wp.site.u_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.v,sol.vv] = deal( Wp.site.v_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.p,sol.pp] = deal( Wp.site.p_init * ones(Wp.mesh.Nx,Wp.mesh.Ny) ); 

    % Initialize the linearized solution variables, if necessary
    if options.Linearversion
        sol.ul = sol.u;
        sol.vl = sol.v;
        sol.pl = sol.p;
        [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
    end;

    % Compute boundary conditions and system matrices B1, B2.
    [sys.B1,sys.B2,sys.bc]  = Compute_B1_B2_bc(Wp);
    sys.pRCM                = []; % Load empty vector
    
    % Compute projection matrices Qsp and Bsp. These are only necessary if
    % the continuity equation is projected away (2015 ACC paper, Boersma).
    if options.Projection
        [sys.Qsp, sys.Bsp]  = Solution_space(sys.B1,sys.B2,sys.bc);
        Wp.Nalpha           = size(sys.Qsp,2);
    end
end

