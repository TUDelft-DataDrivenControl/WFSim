function [ sol_out,sys_out ] = WFSim_timestepping( sol_in,sys_in,Wp,scriptOptions )
    % Import necessary information from sol_in (previous timestep)
    sol_out = struct('k',sol_in.k+1,'u',sol_in.u,'uu',sol_in.uu,...
                     'v',sol_in.v,'vv',sol_in.vv);
    if scriptOptions.exportPressures
        sol_out.p = sol_in.p;
        sol_out.pp = sol_in.pp;
    end;
    
    sys_out = struct('B1',sys_in.B1,'B2',sys_in.B2,'bc',sys_in.bc,...
                     'pRCM',sys_in.pRCM);    
    if scriptOptions.Projection
        sys_out.Qsp = sys_in.Qsp;
        sys_out.Bsp = sys_in.Bsp;
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

        %display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
    end
end

