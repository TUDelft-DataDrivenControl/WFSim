function [ sol_out,sys_out ] = WFSim_timestepping( sol_in, sys_in, Wp, turbInput, modelOptions )
%WFSIM_TIMESTEPPING Propagates the WFSim solution one timestep forward

    % Import necessary information from sol_in (previous timestep)
    sol_out      = sol_in;
    sol_out.k    = sol_in.k + 1;            % Propagate forward in time
    sol_out.time = sol_in.time + Wp.sim.h; % Propagate forward in time
    sol_out.uk   = sol_in.u;
    sol_out.vk   = sol_in.v;
    
    % Write control settings to the solution
    if sol_in.k == 0
        turbInput.dCT_prime = zeros(Wp.turbine.N,1);
    else
        turbInput.dCT_prime = turbInput.CT_prime - sol_in.turbInput.CT_prime;
    end    
    sol_out.turbInput = turbInput;   
    
    % Copy relevant system matrices from previous time
    sys_out = struct('B1',sys_in.B1,'B2',sys_in.B2,'bc',sys_in.bc,'pRCM',sys_in.pRCM); 
    if modelOptions.Projection
        sys_out.Qsp = sys_in.Qsp;
        sys_out.Bsp = sys_in.Bsp;
    end;
    
    % Check if pRCM is properly defined
    if length(sys_out.pRCM) <= 0 && sol_out.k > 1
        error('pRCM not assigned. Please run [sol,sys] = WFSim_timestepping(..) at sol.k = 1 to determine pRCM.'); 
    end
    
    % Load convergence settings
    conv_eps = modelOptions.conv_eps;
    if sol_out.k > 1
        max_it = modelOptions.max_it_dyn;
    else
        max_it = modelOptions.max_it;
    end

    % Initialize default convergence parameters
    it         = 0;
    eps        = 1e19;
    epss       = 1e20;
    
    % Convergence until satisfactory solution has been found
    while ( eps>conv_eps && it<max_it && eps<epss )
        it   = it+1;
        epss = eps;

        % Create system matrices sys.A and sys.b for our nonlinear model,
        % where WFSim basically comes down to: sys.A*sol.x=sys.b.
        [sol_out,sys_out] = Make_Ax_b(Wp,sys_out,sol_out,modelOptions); 
        
        % Compute solution sol.x for sys.A*sol.x = sys.b.
        [sol_out,sys_out] = Computesol(sys_out,sol_out,it,modelOptions); 
        
        % Map solution sol.x to the flow field sol.u, sol.v, sol.p.
        [sol_out,eps] = MapSolution(Wp,sol_out,it,modelOptions); 
    end
end

