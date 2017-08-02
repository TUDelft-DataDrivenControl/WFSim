function [ sol_out,sys_out ] = WFSim_timestepping( sol_in, sys_in, Wp, scriptOptions )

    % Import necessary information from sol_in (previous timestep)
    sol_out      = sol_in;
    sol_out.k    = sol_in.k + 1;            % Propagate forward in time
    sol_out.time = Wp.sim.time(sol_in.k+2); % Propagate forward in time
    sys_out      = struct('B1',sys_in.B1,'B2',sys_in.B2,'bc',sys_in.bc,...
                          'pRCM',sys_in.pRCM);   
                 
    if scriptOptions.Projection
        sys_out.Qsp     = sys_in.Qsp;
        sys_out.Bsp     = sys_in.Bsp;
    end;
    
    % Load convergence settings
    conv_eps    = scriptOptions.conv_eps;
    max_it      = scriptOptions.max_it;
    max_it_dyn  = scriptOptions.max_it_dyn;
    
    % Initialize default convergence parameters
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
 
    while ( eps>conv_eps && it<max_it && eps<epss )
        it   = it+1;
        epss = eps;

        if sol_out.k > 1
            max_it = max_it_dyn;
        end

        [sol_out,sys_out] = Make_Ax_b(  Wp,sys_out,sol_out,   scriptOptions); % Create system matrices
        [sol_out,sys_out] = Computesol( Wp,sys_out,sol_out,it,scriptOptions); % Compute solution
        [sol_out,eps]     = MapSolution(Wp,        sol_out,it,scriptOptions); % Map solution to field
    end
end

