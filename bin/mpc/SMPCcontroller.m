function [CT_prime,phi,mpc] = SMPCcontroller(sol,Wp,NN)

% initialize controller
mpc = MPCinit(sol,Wp,NN);

% solve mpc
if sol.k>=1
    
    xinit         = zeros(mpc.nx*Wp.turbine.N,1);
    xinit(mpc.Mf) = sol.turbine.force;
    xinit(mpc.Mp) = sol.turbine.power;
    xinit(mpc.Mu) = sol.turbine.CT_prime;
    
    yalmip('clear');
    
    cons = [];
    cost = 0;
    
    % define decision variables for the windfarm
    U    = sdpvar(Wp.turbine.N*mpc.Nh,1);
    
    % Ns is number of samples taken into account with the SMPC
    Ns      = 6;
    % mean_v is the measured wind speed
    mean_v  = sol.turbine.Ur;
    % sigma_v is the standard deviation of stochastic wind speed
    sigma_v = .5;
    
    for ll=1:Ns
        
        
        noise     = sigma_v*randn(Wp.turbine.N,mpc.Nh);
        mpc.V     = repmat(mean_v,1,mpc.Nh) + noise;
        
        % build wind farm model
        mpc       = wfmodel(sol,Wp,mpc,2);
        
        % build matrices horizon
        mpc       = matsys(Wp,mpc);
        
        X         = mpc.AA*xinit + mpc.BBt*U ;
        Y         = mpc.CC*X;
        P         = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;
        
        E(:,ll)   = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-sum(P)';
        
        cost = cost + E(:,ll)'*mpc.Q/Ns*E(:,ll);
        
    end
    
    cons = [cons, mpc.um <= U <= mpc.uM];
    
    dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];
    
    cost = cost + dU'*mpc.R*dU;
    
    ops  = sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
    
    optimize(cons,cost,ops)
    
end

%% Assign the decision variables
Yopt          = value(Y);
Uopt          = Yopt(mpc.MU);
temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);

CT_prime      = temp(:,1);              % first action horizon
phi           = zeros(Wp.turbine.N,1);

end