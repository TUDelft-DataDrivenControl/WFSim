function [CT_prime,phi,mpc] = NMPCcontroller(sol,Wp,NN)

% initialize controller
mpc = MPCinit(sol,Wp,NN);

% solve mpc
if sol.k>=1
    
    xinit         = zeros(mpc.nx*Wp.turbine.N,1);
    xinit(mpc.Mf) = sol.turbine.force;
    xinit(mpc.Mp) = sol.turbine.power;
    xinit(mpc.Mu) = sol.turbine.CT_prime;
    
    % nl-1 is number of times the rotor-averaged wind speeds in the horizon
    % will be updated during one sample. If nl=1, the rotor-averaged wind speeds
    % are taken constant in the horizon
    nl = 2;
    
    for ll=1:nl
        
        yalmip('clear');
        
        cons = [];
        
        % define decision variables for the windfarm
        U    = sdpvar(Wp.turbine.N*mpc.Nh,1);
        
        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,ll);
        
        % build matrices horizon
        mpc  = matsys(Wp,mpc);
        
        X    = mpc.AA*xinit + mpc.BBt*U ;
        Y    = mpc.CC*X;
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;
        
        E    = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-sum(P)';
        
        cons = [cons, mpc.um <= U <= mpc.uM];
        
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];
        
        cost = E'*mpc.Q*E + dU'*mpc.R*dU;
        
        ops  = sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
        
        optimize(cons,cost,ops)
        
        Uopt(:,ll) = value(U); %value(Y(mpc.MU))
        Popt       = value(P);
        
        temp       = ( repmat(Popt(:,ll),mpc.Nh,1)./(mpc.cp(1).*Uopt(:,ll))  ).^(1/3);
        mpc.V      = reshape(temp,Wp.turbine.N,mpc.Nh); % rotor-averaged wind speed in the horizons
        
    end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(mpc.MU);
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
end