%
omega = logspace(-3,log10(pi),100)';

for om=1:length(omega)
    for f=1:size(Btl,2);
        FRF(om,:,f)=    full(Qsp)*((Etl*exp(1i*omega(om))-Atl)\Btl(:,f));
    end
end

FRFu   = zeros(size(FRF,1),n2,Wp.Nx-3,size(FRF,3));
FRFu   = reshape(FRF(:,1:size(C1,1),:),size(FRF,1),n2,Wp.Nx-3,size(FRF,3));
FRFv   = zeros(size(FRF,1),n2,Wp.Nx-2,size(FRF,3));
FRFv   = reshape(FRF(:,size(C1,1)+1:size(C1,1)+size(C2,1),:),size(FRF,1),n2,Wp.Nx-2,size(FRF,3));

ypos  = 2;
in    = 1;

figure(3);clf;
for kk=1:size(FRF,1)
    semilogx(omega,20*log10(abs(FRFu(:,ypos,kk,in))))
    grid;xlabel('\omega [rad/s]');ylabel('|\cdot| dB');
    str = sprintf('xpos = %d , xpos llth turbine = %d',kk,Wp.xline(in));
    title(str);
    drawnow;pause;
end