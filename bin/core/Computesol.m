function [sol,sys] = Computesol(sys,input,sol,k,it,options)

Projection      = options.Projection;
exportLinearSol = options.exportLinearSol;

if  Projection
    if (it==1 && k==1)
        % Define initial condition
        sol.alpha = sys.Qsp\([vec(sol.u(3:end-1,2:end-1)');vec(sol.v(2:end-1,3:end-1)')]-sys.Bsp);
    end
    
    Ft                    = sys.At*sol.alpha + sys.Bt*[input.beta;input.phi] + sys.St;
    sol.alpha(sys.pRCM,1) = sys.Et(sys.pRCM,sys.pRCM)\Ft(sys.pRCM);
    sol.x                 = sys.Qsp*sol.alpha + sys.Bsp;
    
    if  exportLinearSol
        if (it==1 && k==1);sol.dalpha = zeros(size(sys.Qsp,2),1);end
        Ftl                     = sys.Atl*sol.dalpha + sys.Btl*[input.dbeta;input.dphi];
        sol.dalpha(sys.pRCM,1)  = sys.Etl(sys.pRCM,sys.pRCM)\Ftl(sys.pRCM);
        sol.dx                  = sys.Qsp*sol.dalpha;
    end
    
else
    sol.x(sys.pRCM,1) = sys.A(sys.pRCM,sys.pRCM)\sys.b(sys.pRCM);
    
    if  exportLinearSol
        if (it==1 && k==1); sol.dx = zeros(size(sys.Al,2),1);end
        bll                 = sys.Al*sol.dx+sys.bl;
        sol.dx(sys.pRCM,1)  = sys.A(sys.pRCM,sys.pRCM)\bll(sys.pRCM);
    end
end