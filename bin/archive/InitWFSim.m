function [Wp,Power,CT,Ueffect,u,v,p,uu,vv,pp,uk,vk,pk,du,dv,dp,input,B1,B2,bc,Qsp,Bsp,solnew] = InitWFSim(Wp)


[Wp,input]   = meshing2(Wp,0,0);         % Create wind farm

% Initial flow fields
[u,uu]      = deal(Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
uk          = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
[v,vv]      = deal(Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
vk          = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
[p,pp]      = deal(Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny)); 
pk          = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
[du,dv,dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));

% Init signals
[Power,CT,Ueffect] = deal(zeros(Wp.turbine.N,Wp.sim.NN)); 

[B1,B2,bc]         = Compute_B1_B2_bc(Wp.mesh,Wp.site.u_Inf);

if nargout>20;
        [Qsp, Bsp]  = Solution_space(B1,B2,bc); % Projection matrices
        solnew      = Qsp\([vec(u(3:end-1,2:end-1)');vec(v(2:end-1,3:end-1)')]...
            -Bsp);                              % Initial condition
end

end

