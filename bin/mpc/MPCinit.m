function mpc = MPCinit(sol,Wp,NN)

% controller parameters
mpc.Nh         = 10;                              % prediction horizon
mpc.Q          = 1e-4*eye(mpc.Nh);                   % weigth on tracking
mpc.R          = 1e9*eye(mpc.Nh*Wp.turbine.N);      % weigth on control signal
mpc.duc        = 1e-1;                            % limitation on du/dt
mpc.um         = .1;                              % minimum CT'
mpc.uM         = 2;                               % maximum CT'

% boleans to extract desired signals from state vector
mpc.MF         = logical(repmat(repmat([1 0 0]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mf         = logical(repmat([1 0 0]',Wp.turbine.N,1));
mpc.MP         = logical(repmat(repmat([0 1 0]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mp         = logical(repmat([0 1 0]',Wp.turbine.N,1));
mpc.MU         = logical(repmat(repmat([0 0 1]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mu         = logical(repmat([0 0 1]',Wp.turbine.N,1));

% wind farm reference
load('bin/mpc/P_reference')
N0                  = round(.1*NN/(2*Wp.sim.h));
mpc.Pref            = zeros(NN+mpc.Nh,1); 
mpc.AGCdata         = AGCdata(:,2); 

mpc.Pgreedy         = 3.438550894184890e+07;                  % 9turb with CT'=2 in steady-state  
mpc.Pref(1:N0)      = .9*mpc.Pgreedy; 
mpc.Pref(N0+1:end)  = .9*mpc.Pgreedy + .2*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);

% controller models
mpc.tau          = 5;                             % time constant filter CT'
mpc.cf           = sol.turbine.cf;
mpc.cp           = sol.turbine.cp;

[num,den]       = tfdata(c2d(tf(1,[mpc.tau 1]),Wp.sim.h),'v');    % filter on force

for kk = 1:Wp.turbine.N    
    mpc.a{kk}     = kron(eye(3),-den(2));
    mpc.bcoef{kk} = num(2)*[-mpc.cf(kk);mpc.cp(kk);1];
    mpc.c{kk}     = eye(3);                                     
end    
mpc.nx = size(mpc.a{1},1);
