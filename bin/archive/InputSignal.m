Phi             = zeros(Wp.N,NN);           % Yaw angles in degrees (-90 < Phi < 90 degrees)
dPhi            = zeros(Wp.N,NN);
beta            = 1/3*ones(Wp.N,NN);        % Scaled axial induction 0 < beta < 1
dbeta           = zeros(Wp.N,NN);


%Phi(1,:)  = 35*[zeros(1,ceil(.1*L)) ones(1,NN-ceil(.1*L))];
%Phi(2,:)  = -5*[zeros(1,ceil(.1*L)) ones(1,NN-ceil(.1*L))];
%Phi(1,:)  = [zeros(1,ceil(.01*L)) 30*ones(1,ceil(.5*NN-ceil(.01*L))) zeros(1,NN-ceil(.5*NN-ceil(.01*L))-ceil(.01*L))];
Phi(1,:)    = [Phi(1,1) Phi(1,2:end)+0*ones(1,NN-1)];
Phi(2,:)    = [Phi(2,1) Phi(2,2:end)-0*ones(1,NN-1)];
beta(1,:)   = [beta(1,1) beta(1,2:end)+.05*ones(1,NN-1)];
beta(2,:)   = [beta(2,1) beta(2,2:end)-.05*ones(1,NN-1)];

%tau1 = .05;
%Phi  = lsim(ss(-1/tau1,1/tau1,1,0)*eye(Wp.N),Phi,time,Phi(:,1))';

%for kk=1:NN;if sin(2*pi*2/L*time(kk))>=0;beta(:,kk) = .5;else beta(:,kk) = .2;end;end;

%tau2 = 10;
%beta = lsim(ss(-1/tau2,1/tau2,1,0)*eye(Wp.N),beta,time,beta(:,1))';

dbeta(:,:) = [diff(beta')' diff(beta(:,end-1:end)')'];
dPhi(:,:)  = [diff(Phi')' diff(Phi(:,end-1:end)')'];

