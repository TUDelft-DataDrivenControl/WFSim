function [sol,eps] = MapSolution(Wp,sol,it,options)
%MAPSOLUTION This function converts 'sol.x' to real flow fields.

    % Import variables
    k  = sol.k;
    Nx = Wp.mesh.Nx;
    Ny = Wp.mesh.Ny;
    exportPressures   = options.exportPressures;
    exportLinearSol   = options.exportLinearSol;

    % Project sol.x back to the flow fields, excluding the boundary conditions
    sol.uu(3:end-1,2:end-1)     = reshape(sol.x(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
    sol.vv(2:end-1,3:end-1)     = reshape(sol.x((Nx-3)*(Ny-2)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)),Ny-3,Nx-2)';
    %pp(2:end-1,2:end-1)     = reshape([sol.x((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0],Ny-2,Nx-2)';
    if exportPressures
        sol.pp(2:end-1,2:end-1) = reshape([sol.x((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0;0],Ny-2,Nx-2)';
        sol.pp(isinf(sol.pp))   = 0;
    end

    if exportLinearSol
        sol.du(3:end-1,2:end-1) = reshape(sol.dx(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
        sol.dv(2:end-1,3:end-1) = reshape(sol.dx((Nx-3)*(Ny-2)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)),Ny-3,Nx-2)';
        sol.ul(3:end-1,2:end-1) = sol.ul(3:end-1,2:end-1)+sol.du(3:end-1,2:end-1);
        sol.vl(2:end-1,3:end-1) = sol.vl(2:end-1,3:end-1)+sol.dv(2:end-1,3:end-1);
        if exportPressures
            %sol.dp(2:end-1,2:end-1)     = reshape([sol.dx((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0],Ny-2,Nx-2)';
            sol.dp(2:end-1,2:end-1)     = reshape([sol.dx((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0;0],Ny-2,Nx-2)';
            sol.dp(isinf(sol.dp))       = 0;

            %sol.pl(2:end-1,2:end-1)     = sol.pl(2:end-1,2:end-1)+sol.dp(2:end-1,2:end-1);
            sol.pl(2:end-1,2:end-1)     = reshape([sol.dx((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0;0],Ny-2,Nx-2)';
            sol.pl(isinf(sol.pl))       = 0;
        end
    end

    % Check if solution has converged
    Normv = norm(vec(sol.v(2:end-1,3:end-1)-sol.vv(2:end-1,3:end-1)));
    Normu = norm(vec(sol.u(3:end-1,2:end-1)-sol.uu(3:end-1,2:end-1)));
    eps   = sqrt((Normv+Normu))/((Ny-2)*(Nx-2))/2;

    if k==1; alpha = min(1-.9^it,1); else; alpha=1; end;
    sol.u(3:end-1,2:end-1)  = (1-alpha)*sol.u(3:end-1,2:end-1)+alpha*sol.uu(3:end-1,2:end-1);
    sol.v(2:end-1,3:end-1)  = (1-alpha)*sol.v(2:end-1,3:end-1)+alpha*sol.vv(2:end-1,3:end-1);
    sol.p(2:end-1,2:end-1)  = (1-alpha)*sol.p(2:end-1,2:end-1)+alpha*sol.pp(2:end-1,2:end-1);

    % Update velocities for next iteration and boundary conditions
    [sol.u,sol.v,sol.p] = Updateboundaries(Nx,Ny,sol.u,sol.v,sol.p);
    if exportLinearSol
        [sol.ul,sol.vl,sol.pl] = Updateboundaries(Nx,Ny,sol.ul,sol.vl,sol.pl);
        [sol.du,sol.dv,sol.dp] = Updateboundaries(Nx,Ny,sol.du,sol.dv,sol.dp);
    end

    if options.printConvergence
        disp(['k ',num2str(sol.k,'%-1000.f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv,'%10.2e'), ', Nu=', num2str(Normu,'%10.2e'), ', TN=',num2str(eps,'%10.2e'),'']) ;
    end;
end