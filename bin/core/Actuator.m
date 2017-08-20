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
F               = Wp.turbine.forcescale; % http://www.nrel.gov/docs/fy05osti/36834.pdf
input           = Wp.turbine.input(sol.k);
Projection      = options.Projection;
Linearversion   = options.Linearversion;
Derivatives     = options.Derivatives;

%%
Ar              = pi*(0.5*Drotor)^2;
scale           = 3;                                  % To scale the force in the y-direction

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
    % Clear for multiple turbine case
    tempx = sparse(Nx-3,Ny-2);
    tempy = sparse(Nx-2,Ny-3);
    
    x  = xline(kk,:);  % Turbine x-pos in field
    y  = yline{kk};    % Turbine y-pos in field
    yv = ylinev{kk};   % Corrected turbine y-pos in field
    
    vv            = 0.5*diff(sol.v(x,yv))+sol.v(x,yv(1:end-1)); % Bart: this can be fixed!
    uu            = sol.u(x,y);
    U{kk}         = sqrt(uu.^2+vv.^2);
    phi{kk}       = atan(vv./uu);
    Ue{kk}        = cos(phi{kk}+input.phi(kk)/180*pi).*U{kk};
    meanUe{kk}    = mean(Ue{kk});
    dUedPhi{kk}   = -1/180*pi*sin(phi{kk} + input.phi(kk)/180*pi).*U{kk} ;
    CT(kk)        = input.CT_prime(kk); % Import CT_prime from inputData
       
 %% Thrust force      
    % With the following, we only take middle velocity component and
    % apply this on the rotor cells 
    m = size(Ue{kk},2);
    if rem(m,2)
       ind  = ceil(m/2); 
       temp = Ue{kk};
       Ur   = repmat(temp(ind),1,m);
    else
       ind  = [m/2 m/2+1]; 
       temp = Ue{kk};
       Ur   = repmat(mean(temp(ind)),1,m);
    end    
    Fthrust         = F*1/2*Rho*Ue{kk}.^2*CT(kk); % Using CT_prime
    Fx              = Fthrust.*cos(input.phi(kk)*pi/180);
    Fy              = Fthrust.*sin(input.phi(kk)*pi/180);
    
    %%
    pp = 1.88; % Loss factor for yawing a turbine
    Power(kk) = powerscale*mean(.5*Rho*Ar*Ue{kk}.^3*CT(kk)*cos(input.phi(kk)*pi/180)^pp);    

    %% Input to Ax=b
    Sm.x(x-2,y-1)           = -Fx'.*dyy2(1,y)';                                                                  % Input x-mom nonlinear                           
    Sm.y(x-1,y(2:end)-2)    = scale*Fy(2:end)'.*dyy2(1,y(2:end))';                                               % Input y-mom nonlinear
    
    if Linearversion
        tempdx = sparse(Nx-3,Ny-2);
        tempdy = sparse(Nx-2,Ny-3);
        
        dFthrustdbeta   = 1/2*Rho*Ue{kk}.^2*(CT(kk)*2*(input.beta(kk)+1)+dCTdbeta(kk)*(input.beta(kk)+1).^2);
        dFxdbeta        = dFthrustdbeta.*cos(input.phi(kk)*pi/180);
        dFydbeta        = dFthrustdbeta.*sin(input.phi(kk)*pi/180);
        
        %% Gradient of the thrust with respect to yaw angle
        dFthrustdPhi    = Rho*CT(kk)*Ue{kk}.*(input.beta(kk)+1).^2.*dUedPhi{kk};
        dFxdPhi         = dFthrustdPhi.*cos(input.phi(kk)*pi/180+0*phi{kk}) - pi/180*Fthrust.*sin(input.phi(kk)*pi/180);
        dFydPhi         = dFthrustdPhi.*sin(input.phi(kk)*pi/180+0*phi{kk}) + pi/180*Fthrust.*cos(input.phi(kk)*pi/180);
        
        Sm.dx(x-2,y-1)          = -(dFxdbeta'*input.dbeta(kk) + dFxdPhi'*input.dphi(kk)).*dyy2(1,y)';
        Sm.dy(x-1,y(2:end)-2)   = scale*(dFydbeta(2:end)'*input.dbeta(kk) + dFydPhi(2:end)'*input.dphi(kk)).*dyy2(1,y(2:end))';  % Input y-mom linear
        
        tempdx(x-2,y-1) = -dFxdbeta'.*dyy2(1,y)';
        Sm.dxx(:,kk)    =  vec(tempdx');           % Input matrix (beta) x-mom linear
        tempdx(x-2,y-1) = -dFxdPhi'.*dyy2(1,y)';
        Sm.dxx(:,N+kk)  =  vec(tempdx');           % Input matrix (yaw) x-mom linear
        
        tempdy(x-1,y(2:end)-2)  = dFydbeta(2:end)'.*dyy2(1,y(2:end))';
        Sm.dyy(:,kk)            = scale*vec(tempdy');                   % Input (beta) y-mom linear qlpv
        tempdy(x-1,y(2:end)-2)  = dFydPhi(2:end)'.*dyy2(1,y(2:end))';
        Sm.dyy(:,N+kk)          = scale*vec(tempdy');                   % Input (yaw) y-mom linear qlpv
    end;
    
    if Projection
        %% Input to qLPV
        tempx(x-2,y-1)          = -Fx'.*dyy2(1,y)';
        Sm.xx(:,kk)             = vec(tempx')/input.beta(kk);
        Sm.xx(:,N+kk)           = vec(tempx');
        if input.phi(kk)~=0
            Sm.xx(:,kk)         = Sm.xx(:,kk)/2;
            Sm.xx(:,N+kk)       = Sm.xx(:,N+kk)/(2*input.phi(kk));         % Input x-mom nonlinear qlpv
        end
        
        tempy(x-1,y(2:end)-2)   = Fy(2:end)'.*dyy2(1,y(2:end))';
        Sm.yy(:,kk)             = scale*vec(tempy')/input.beta(kk);
        Sm.yy(:,N+kk)           = scale*vec(tempy');
        if input.phi(kk)~=0
            Sm.yy(:,kk)         = Sm.yy(:,kk)/2;
            Sm.yy(:,N+kk)       = Sm.yy(:,N+kk)/(2*input.phi(kk));         % Input y-mom nonlinear qlpv
        end
        
    end;
    
    if Linearversion
        dFu           = 4*input.beta(kk)*Rho*uu*cos(input.phi(kk)*pi/180);
        dFv           = 4*input.beta(kk)*Rho*sol.v(x,y)*sin(input.phi(kk)*pi/180);      % SB: why here original v and not vv?
        Smdu(x-2,y-1) = -dFu'.*dyy2(1,y)';
        Smdv(x-1,y-2) =  dFv'.*dyy2(1,y)';
        
        dSm.dx = blkdiag(diag(vec(Smdu')'),diag(vec(Smdv')'),sparse((Ny-2)*(Nx-2),(Ny-2)*(Nx-2)));
    end;
    
    if Derivatives
        
        dvvdv         = 1/2*diag(sign(ylinev{kk}))+1/2*diag(sign(ylinev{kk}(1:end-1)),1); dvvdv(end,:)=[];
        dUdv(kk,:)    =  vv./U{kk};
        dUdu(kk,:)    =  uu./U{kk};
        
        dUedv(kk,:)=-sin(atan(vv./(uu)) + input.phi(kk)/180*pi)./((uu).*...
            (vv.^2./(uu).^2 + 1)).*U{kk}+cos(atan(vv./(uu))+input.phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
        dUedv(kk,:)=-uu.*sin(atan(vv./uu) + input.phi(kk)/180*pi)./(vv.^2+uu.^2).*U{kk}+cos(atan(vv./(uu))+...
            input.phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
        dUedu(kk,:)=(vv.*sin(atan(vv./(uu)) + input.phi(kk)/180*pi))./((vv.^2 + (uu).^2)).*...
            U{kk}+cos(atan(vv./(uu))+input.phi(kk)/180*pi).*dUdu(kk,:); %% diff Ue wrt u
        
        % Derivatives Thrust with respect to the state
        dFxdU = Rho*CT(kk)*Ue{kk}*((input.beta(kk)+1)).^2*cos((mean(0*phi{kk})+input.phi(kk))*pi/180);
        dFydU = Rho*CT(kk)*Ue{kk}*((input.beta(kk)+1)).^2*sin((mean(0*phi{kk})+input.phi(kk))*pi/180);
        dSm.xdu((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-3)*(Ny-2)+yline{kk}-1) = diag((-dFxdU.*dUedu(kk,:).*dyy2(1,yline{kk})));
        %dSm.xdv((xline(kk,:)-3)*(Ny-2)+ylinev{kk}-1,(xline(kk,:)-2)*(Ny-3)+ylinev{kk}-2) = diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})*dvvdv));
        dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-2) = 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
        dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1) = dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1)+ 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
        dSm.ydu((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-3)*(Ny-2)+yline{kk}(2:end)-1) = diag(scale*(dFydU(2:end).*dUedu(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        %dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}-2) = diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))*dvvdv(2:end,2:end)));
        dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2) = 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1) =   dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1)+ 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
        
        
        dSm.dbetadPower_in(kk) = 1/(1/F*powerscale*2*Rho*pi*(0.5*Drotor)^2*mean(Ue{kk}).^3);
        dSm.dJdPower_in(kk)    = -1/F*powerscale*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3*dSm.dbetadPower_in(kk);
        
        
        dSm.dJdPower_in(kk)   = -1/F*powerscale*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3*dSm.dbetadPower_in(kk);
%         if nargin==7
%             JXXu(Wp.xline(kk,:)-2,Wp.yline{kk}-1) =  ((-1/F* Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue{kk}).^2)./(length(Wp.yline{kk}))+ dSm.dJdbeta(kk).*dSm.dbetadUe{kk}).*dUedu{kk};%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
%             JXXv(Wp.xline(kk,:)-1,Wp.ylinev{kk}-2) = ((-1/F*Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue{kk}).^2)./((length(Wp.yline{kk})))+ dSm.dJdbeta(kk).*dSm.dbetadUe{kk}).*dUedv{kk}*dvvdv;%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
%             
%         else
%             JXXu(Wp.xline(kk,:)-2,Wp.yline{kk}-1) =  -1/F* Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue{kk}).^2.*dUedu{kk}./(length(Wp.yline{kk}));%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
%             JXXv(Wp.xline(kk,:)-1,Wp.ylinev{kk}-2) = -1/F*Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue{kk}).^2.*dUedv{kk}./((length(Wp.yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
%         end
      
        % Derivative cost function with respect to state and beta
        JXXu(xline(kk,:)-2,yline{kk}-1) =  -1/F*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue{kk}).^2.*dUedu(kk,:)./(length(yline{kk}));%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
        JXXv(xline(kk,:)-1,ylinev{kk}-2) = -1/F*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue{kk}).^2.*dUedv(kk,:)./((length(yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
        
        % two lines above should be checked!!
        dSm.dJdx        = [vec(JXXu');vec(JXXv');vec(JXXp')];
        dSm.dJdbeta(kk) = -1/F*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue{kk}).^3;
        dSm.dJdPhi(kk)  = -1/F*6*Rho*pi*(0.5*Drotor)^2*input.beta(kk)*mean(Ue{kk}).^2.*mean(dUedPhi{kk});
        
        %
        dSmxdbetakk = sparse(Nx-3,Ny-2);
        dSmydbetakk = sparse(Nx-2,Ny-3);
        dSmxdphikk  = sparse(Nx-3,Ny-2);
        dSmydphikk  = sparse(Nx-2,Ny-3);
        
        dSmxdbetakk(xline(kk,:)-2,yline{kk}-1)        =  -(dFxdbeta.*dyy2(1,yline{kk}))';  %JW changed definition of force and location (need correction for yaw
        dSmydbetakk(xline(kk,:)-1,yline{kk}(2:end)-2) = scale*(dFydbeta(2:end).*dyy2(1,yline{kk}(2:end)))'; %% JW made this 0
        dSm.dbeta(:,kk)                               = [vec(dSmxdbetakk');vec(dSmydbetakk')];
        
        dSmxdphikk(xline(kk,:)-2,yline{kk}-1)         = -dFxdPhi'.*dyy2(1,yline{kk})';
        dSmydphikk(xline(kk,:)-1,yline{kk}(2:end)-2)  = scale*dFydPhi(2:end)'.*dyy2(1,yline{kk}(2:end))';
        dSm.dphi(:,kk)                                = [vec(dSmxdphikk');vec(dSmydphikk')];
        
        dSm.dPower_in(:,kk)                           = [vec(dSmxdbetakk');vec(dSmydbetakk')].*dSm.dbetadPower_in(kk);

    end
end

%% Write to outputs
sol.turbine.power    = Power;
sol.turbine.CT_prime = CT;

output.Sm  = Sm;
if (Derivatives>0 || Linearversion>0)
    output.dSm = dSm;
end