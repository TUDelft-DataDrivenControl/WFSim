% After running WFSim with the MPC, you can run this function to analyse
% the results

for kk=1:sol.k
    
    time(kk)       = sol_array(kk).k;
    Power(:,kk)    = sol_array(kk).turbine.power;
    CT_prime(:,kk) = sol_array(kk).turbine.CT_prime;
    e(kk)          = mpc.Pref(kk)-sum(Power(:,kk)); 
end

figure(2);clf
subplot(2,1,1)
stairs(time,sum(Power)/1e6,'linewidth',2);hold on;grid
stairs(time,mpc.Pref(1:sol.k)/1e6,'r--','linewidth',2);
ylim([.8*min(mpc.Pref(1:sol.k))/1e6 1.2*max(mpc.Pref(1:sol.k))/1e6])
xlabel('$k$','interpreter','latex')
ylabel('$\sum P_i$ [MW]','interpreter','latex');
title('Wind farm power (blue) and reference (red)','interpreter','latex');
subplot(2,1,2)
for kk=1:Wp.turbine.N
stairs(time,CT_prime(kk,:),'linewidth',2);hold on;
end
grid;ylim([0 mpc.uM])
xlabel('$k$','interpreter','latex')
ylabel('$C_T^\prime$','interpreter','latex');

