for kk = 1:1179
    load(strcat('Run_',num2str(kk)));    
    D(kk,:) = [VAF_1';VAF_2';VAF_3';RMSE_power';VAF_power'];
end

Dmax = max(D);