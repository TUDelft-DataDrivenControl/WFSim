function [B1,B2,Bm1,Bm2,bc] = Compute_B1_B2_bc(Wp,u_Inf,Rho)
Bm1                  = Rho*(spdiags(-ones((Wp.Nx-2)*(Wp.Ny-2),1).*vec(Wp.dyy2(2:end-1,2:end-1)'),0,(Wp.Nx-3)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2))+...
    spdiags(ones((Wp.Nx-2)*(Wp.Ny-2),1).*vec(Wp.dyy2(2:end-1,2:end-1)'),Wp.Ny-2,(Wp.Nx-3)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2)));
Bm2                  = Rho*(spdiags(-ones((Wp.Nx-2)*(Wp.Ny-2),1).*vec(Wp.dxx2(2:end-1,2:end-1)'),0,(Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2))+...
    spdiags(ones((Wp.Nx-2)*(Wp.Ny-2),1).*vec(Wp.dxx2(2:end-1,2:end-1)'),1,(Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2)));
Bm2(Wp.Ny-2:Wp.Ny-2:end,:) = [];

B1 = Bm1'; B2 = Bm2';

bc                  = zeros((Wp.Ny-2)*(Wp.Nx-2),1);
bc(1:Wp.Ny-2)       = -Rho*u_Inf*Wp.dyy2(1,2:end-1)';

B1((Wp.Ny-2)*(Wp.Nx-3)+1:(Wp.Ny-2)*(Wp.Nx-2),:) = 0; % u_{Wp.Nx,J}=u_{Wp.Nx-1,J}


%% Lines 23-30 slow down WFSim/WFObs by changing the system sparsity.
% Thus, we can no longer use the RCM algorithm efficiently... Compare
% spy(Et) with spy(Et(pRCM,pRCM)) to visualize this effect. However,
% simulations show the neccesity of these lines for good estimations.

for kk=1:Wp.Nx-2
    B2(kk*(Wp.Ny-2),:)= 0; % v_{I,Wp.Ny}=v_{I,Wp.Ny-1}
end

% Next lines should be there but do not influence the solution
for kk=0:Wp.Nx-3
    B2(kk*(Wp.Ny-2)+1,:)= 0; % v_{I,3}=v_{I,2} for I=2,3,...,Wp.Nx
end

%%

B1 = B1';
B2 = B2';