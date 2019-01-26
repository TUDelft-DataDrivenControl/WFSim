function mpc = wfmodel(sol,Wp,mpc,ll)

% build B matrix of the turbine models for current times step
for kk = 1:Wp.turbine.N
    
    if ll==1
        v           = repmat(sol.turbine.Ur(kk),mpc.Nh,1);
    else
        v           = mpc.V(kk,:);
    end
    mpc.b{kk}   = [mpc.bcoef{kk}(1)*v(1)^2 mpc.bcoef{kk}(2)*v(1)^3 mpc.bcoef{kk}(3)]';
    
    for nn = 1:mpc.Nh
        mpc.bt(:,kk,nn) = [mpc.bcoef{kk}(1)*v(nn)^2 mpc.bcoef{kk}(2)*v(nn)^3 mpc.bcoef{kk}(3)]';
    end
    
end

% build wind farm model
mpc.A = blkdiag(mpc.a{:});
mpc.B = blkdiag(mpc.b{:});
for kk = 1:mpc.Nh
    for nn = 1:Wp.turbine.N
        b{nn}     = mpc.bt(:,nn,kk);
    end
    mpc.Bt(:,:,kk) = blkdiag(b{:});
end
mpc.C = blkdiag(mpc.c{:});

end

