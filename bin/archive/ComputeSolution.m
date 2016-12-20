if  Projection==1;

    Ct                      = blkdiag(spdiags(ccx,0,Wp.nu,Wp.nu),spdiags(ccy,0,Wp.nv,Wp.nv));
    Et                      = Qsp'*blkdiag(Ax,Ay)*Qsp;
    
    if (it==1 && k==1);pRCM = symrcm(Et);end
    
    At                      = Qsp'*Ct*Qsp;
    St                      = Qsp'*[bx;by] - Qsp'*blkdiag(Ax,Ay)*Bsp + Qsp'*Ct*Bsp;
    Bt                      = Qsp'*[Sm.xx;Sm.yy];
    Ft                      = At*alpha + Bt*beta(:,k) + St;
    
    alpha(pRCM,1)           = Et(pRCM,pRCM)\Ft(pRCM);
    sol                     = Qsp*alpha + Bsp;
    uu(3:end-1,2:end-1)     = reshape(sol(1:Wp.nu),Wp.Ny-2,Wp.Nx-3)';
    vv(2:end-1,3:end-1)     = reshape(sol(Wp.nu+1:Wp.nu+Wp.nv),Wp.Ny-3,Wp.Nx-2)';
else
    
    % Remove zero columns and rows from A and b
    A(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    b(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    A(:,end) = [];A(end,:) = [];b(end)=[];
    
    A(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = []; if (it==1 && k==1);pRCM = symrcm(A);end
    sol(pRCM,1)             = A(pRCM,pRCM)\b(pRCM);
    %sol                     = A\b;
    
    uu(3:end-1,2:end-1)     = reshape(sol(1:Wp.nu),Wp.Ny-2,Wp.Nx-3)';
    vv(2:end-1,3:end-1)     = reshape(sol(Wp.nu+1:Wp.nu+Wp.nv),Wp.Ny-3,Wp.Nx-2)';
    %pp(2:end-1,2:end-1)     = reshape([sol(Wp.nu+Wp.nv+1:end);0],Wp.Ny-2,Wp.Nx-2)';
    pp(2:end-1,2:end-1)     = reshape([sol(Wp.nu+Wp.nv+1:end);0;0],Wp.Ny-2,Wp.Nx-2)';
    pp(isinf(pp))           = 0;
end