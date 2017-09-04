function [output,sol] = Actuator(Wp,sol,options)
% Import variables
Nx              = Wp.mesh.Nx;
Ny              = Wp.mesh.Ny;
dyy2            = Wp.mesh.dyy2;
xline           = Wp.mesh.xline;
yline           = Wp.mesh.yline;
ylinev          = Wp.mesh.ylinev;
Rho             = Wp.site.Rho;
Drotor          = Wp.turbine.Drotor;
powerscale      = Wp.turbine.powerscale;
N               = Wp.turbine.N;
F               = Wp.turbine.forcescale; 
input           = Wp.turbine.input(sol.k);
Projection      = options.Projection;
Linearversion   = options.Linearversion;
Derivatives     = options.Derivatives;

%%
Ar              = pi*(0.5*Drotor)^2;

[Sm.x,Sm.dx]    = deal(sparse(Nx-3,Ny-2));            % Input x-mom nonlinear and linear
[Sm.y,Sm.dy]    = deal(sparse(Nx-2,Ny-3));            % Input y-mom nonlinear and linear
[Sm.xx,Sm.dxx]  = deal(sparse((Nx-3)*(Ny-2),2*N));    % Input x-mom nonlinear and linear qlpv
[Sm.yy,Sm.dyy]  = deal(sparse((Nx-2)*(Ny-3),2*N));    % Input y-mom nonlinear and linear qlpv


if Derivatives > 0
    dSm.xdu     = sparse((Nx-3)*(Ny-2),(Nx-3)*(Ny-2));
    dSm.xdv     = sparse((Nx-3)*(Ny-2),(Nx-2)*(Ny-3));
    dSm.ydu     = sparse((Nx-2)*(Ny-3),(Nx-3)*(Ny-2));
    dSm.ydv     = sparse((Nx-2)*(Ny-3),(Nx-2)*(Ny-3));
    JXXu        = sparse(Nx-3,Ny-2);
    JXXv        = sparse(Nx-2,Ny-3);
    JXXp        = sparse(Nx-2,Ny-2);
end

if Linearversion
    Smdu            = sparse(Nx-3,Ny-2);
    Smdv            = sparse(Nx-2,Ny-3);
end;

for kk=1:N
   
    x  = xline(kk,:);  % Turbine x-pos in field
    y  = yline{kk};    % Turbine y-pos in field
    yv = ylinev{kk};   % Corrected turbine y-pos in field
    
    vv            = 0.5*diff(sol.v(x,yv))+sol.v(x,yv(1:end-1)); 
    uu            = sol.u(x,y);
    U{kk}         = sqrt(uu.^2+vv.^2);
    phi{kk}       = atan(sol.v(1,1)/sol.u(1,1));
    Ue{kk}        = cos(input.phi(kk)/180*pi).*U{kk};
    meanUe{kk}    = mean(Ue{kk});  
    CT(kk)        = input.CT_prime(kk); % Import CT_prime from inputData
    dCT(kk)       = input.dCT_prime(kk); 
    Phi(kk)       = input.phi(kk);
    
    %% Thrust force       
    Fthrust         = F*1/2*Rho*Ue{kk}.^2*CT(kk); % Using CT_prime
    Fx              = Fthrust.*cos(phi{kk}+Phi(kk)*pi/180);
    Fy              = Fthrust.*sin(phi{kk}+Phi(kk)*pi/180);
    
    %% Power
    Power(kk)   = powerscale*.5*Rho*Ar*CT(kk)*mean(Ue{kk}.^3);    
    
    %% Input to Ax=b
    Sm.x(x-2,y-1)           = -Fx'.*dyy2(1,y)';               % Input x-mom nonlinear                           
    Sm.y(x-1,y(2:end)-2)    = Fy(2:end)'.*dyy2(1,y(2:end))';  % Input y-mom nonlinear
    
    % Matrices for linear version
    if Linearversion
              
        dFthrustdCT             = F*1/2*Rho*Ue{kk}.^2;
        dFxdCT                  = dFthrustdCT.*cos(phi{kk}+Phi(kk)*pi/180);
        dFydCT                  = dFthrustdCT.*sin(phi{kk}+Phi(kk)*pi/180);
         
        Sm.dx(x-2,y-1)          = -dFxdCT'*dCT(kk).*dyy2(1,y)';
        Sm.dy(x-1,y(2:end)-2)   = dFydCT(2:end)'*dCT(kk).*dyy2(1,y(2:end))';  
               
        dFdu                    = F*Rho*cos(Phi(kk)*pi/180)^2*CT(kk)*uu;
        dFdv                    = F*Rho*cos(Phi(kk)*pi/180)^2*CT(kk)*vv;     
        Smdu(x-2,y-1)           = -dFdu'.*dyy2(1,y)';
        Smdv(x-1,y-2)           =  dFdv'.*dyy2(1,y)';
        
        dSm.dx                  = blkdiag(diag(vec(Smdu')'),diag(vec(Smdv')'),sparse((Ny-2)*(Nx-2),(Ny-2)*(Nx-2)));
        
        % following for projection
        tempdx                  = sparse(Nx-3,Ny-2);
        tempdy                  = sparse(Nx-2,Ny-3);
        tempdx(x-2,y-1)         = -dFxdCT'.*dyy2(1,y)';
        Sm.dxx(:,kk)            =  vec(tempdx');                  % Input matrix (beta) x-mom linear       
        tempdy(x-1,y(2:end)-2)  = dFydCT(2:end)'.*dyy2(1,y(2:end))';
        Sm.dyy(:,kk)            = vec(tempdy');                   % Input (beta) y-mom linear qlpv
    end;
    
    if Projection
        %% Input to qLPV
        % Clear for multiple turbine case
        tempx                   = sparse(Nx-3,Ny-2);
        tempy                   = sparse(Nx-2,Ny-3);
        tempx(x-2,y-1)          = -Fx'.*dyy2(1,y)';
        Sm.xx(:,kk)             = vec(tempx')/CT(kk);
        Sm.xx(:,N+kk)           = vec(tempx');
        if input.phi(kk)~=0
            Sm.xx(:,kk)         = Sm.xx(:,kk)/2;
            Sm.xx(:,N+kk)       = Sm.xx(:,N+kk)/(2*Phi(kk));         % Input x-mom nonlinear qlpv
        end
        
        tempy(x-1,y(2:end)-2)   = Fy(2:end)'.*dyy2(1,y(2:end))';
        Sm.yy(:,kk)             = vec(tempy')/CT(kk);
        Sm.yy(:,N+kk)           = vec(tempy');
        if input.phi(kk)~=0
            Sm.yy(:,kk)         = Sm.yy(:,kk)/2;
            Sm.yy(:,N+kk)       = Sm.yy(:,N+kk)/(2*Phi(kk));         % Input y-mom nonlinear qlpv
        end
        
    end;
      
    % Did not look at derivatives yet
    if Derivatives     
        dvvdv         = 1/2*diag(sign(ylinev{kk}))+1/2*diag(sign(ylinev{kk}(1:end-1)),1); dvvdv(end,:)=[];
        dUdv(kk,:)    =  vv./U{kk};
        dUdu(kk,:)    =  uu./U{kk};
               
        dUedv(kk,:)= cos(input.phi(kk)*pi/180)*v./U{kk};
        dUedu(kk,:)= cos(input.phi(kk)*pi/180)*u./U{kk};
        
        % Derivatives Thrust with respect to the state (SB: are these
        % correct? Not wrt u and v?
        dFxdU = Rho*CT(kk)*Ue{kk}*cos(input.phi(kk)*pi/180);
        dFydU = Rho*CT(kk)*Ue{kk}*sin(input.phi(kk)*pi/180);
              
        dSm.xdu((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-3)*(Ny-2)+yline{kk}-1) = diag((-dFxdU.*dUedu(kk,:).*dyy2(1,yline{kk})));
        dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-2) = 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
        dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1) = dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1)+ 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
        dSm.ydu((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-3)*(Ny-2)+yline{kk}(2:end)-1) = diag((dFydU(2:end).*dUedu(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2) = 1/2*diag((dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1) =   dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1)+ 1/2*diag((dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        
        
        dSm.dbetadPower_in(kk) = 1/(1/cf*cp*2*Rho*pi*(0.5*Drotor)^2*mean(Ue{kk}).^3);
        dSm.dJdPower_in(kk)    = -1/cf*cp*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3*dSm.dbetadPower_in(kk);
        
        
        dSm.dJdPower_in(kk)   = -1/cf*cp*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3*dSm.dbetadPower_in(kk);
      
        % Derivative cost function with respect to state and beta
        JXXu(xline(kk,:)-2,yline{kk}-1) =  -1/cf*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue{kk}).^2.*dUedu(kk,:)./(length(yline{kk}));%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
        JXXv(xline(kk,:)-1,ylinev{kk}-2) = -1/cf*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue{kk}).^2.*dUedv(kk,:)./((length(yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
        
        % two lines above should be checked!!
        dSm.dJdx        = [vec(JXXu');vec(JXXv');vec(JXXp')];
        dSm.dJdbeta(kk) = -1/cf*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3;
        dSm.dJdPhi(kk)  = -1/cf*6*Rho*pi*(0.5*Drotor)^2*input.beta(kk)*mean(Ue{kk}).^2.*mean(dUedPhi{kk});
        
        %
        dSmxdbetakk = sparse(Nx-3,Ny-2);
        dSmydbetakk = sparse(Nx-2,Ny-3);
        dSmxdphikk  = sparse(Nx-3,Ny-2);
        dSmydphikk  = sparse(Nx-2,Ny-3);
        
        dSmxdbetakk(xline(kk,:)-2,yline{kk}-1)        =  -(dFxdCT.*dyy2(1,yline{kk}))';  %JW changed definition of force and location (need correction for yaw
        dSmydbetakk(xline(kk,:)-1,yline{kk}(2:end)-2) = (dFydCT(2:end).*dyy2(1,yline{kk}(2:end)))'; %% JW made this 0
        dSm.dbeta(:,kk)                               = [vec(dSmxdbetakk');vec(dSmydbetakk')];
        
        dSmxdphikk(xline(kk,:)-2,yline{kk}-1)         = -dFxdPhi'.*dyy2(1,yline{kk})';
        dSmydphikk(xline(kk,:)-1,yline{kk}(2:end)-2)  = dFydPhi(2:end)'.*dyy2(1,yline{kk}(2:end))';
        dSm.dphi(:,kk)                                = [vec(dSmxdphikk');vec(dSmydphikk')];
        
        dSm.dPower_in(:,kk)                           = [vec(dSmxdbetakk');vec(dSmydbetakk')].*dSm.dbetadPower_in(kk);

    end
end

%% Write to outputs
sol.turbine.power(:,1)    = Power;
sol.turbine.CT_prime(:,1) = CT;

output.Sm  = Sm;
if (Derivatives>0 || Linearversion>0)
    output.dSm = dSm;
end