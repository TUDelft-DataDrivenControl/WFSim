function [sol,sys] = Computesol(sys,sol,it,options)
% COMPUTESOL  This function solves A*x=b for x to obtain the flow fields.

    % Import variables
    k               = sol.k;
    input           = sol.turbInput;
    Projection      = options.Projection;
    exportLinearSol = options.exportLinearSol;

    % Find the solution
    if  Projection
        % Solve by projecting away the continuity equation
        if (it==1 && k==1)
            % Define initial condition
            sol.alpha = sys.Qsp\([vec(sol.u(3:end-1,2:end-1)');vec(sol.v(2:end-1,3:end-1)')]-sys.Bsp);
        end
        beta                  = 1/4*input.CT_prime;
        Ft                    = sys.At*sol.alpha + sys.Bt*[beta;input.phi] + sys.St;
        sol.alpha(sys.pRCM,1) = sys.Et(sys.pRCM,sys.pRCM)\Ft(sys.pRCM);
        sol.x                 = sys.Qsp*sol.alpha + sys.Bsp;

        if  exportLinearSol
            if (it==1 && k==1);sol.dalpha = zeros(size(sys.Qsp,2),1);end
            dbeta                   = 1/4*input.dCT_prime;
            Ftl                     = sys.Atl*sol.dalpha + sys.Btl*[dbeta;input.dphi];
            sol.dalpha(sys.pRCM,1)  = sys.Etl(sys.pRCM,sys.pRCM)\Ftl(sys.pRCM);
            sol.dx                  = sys.Qsp*sol.dalpha;
        end

    else
        % Otherwise, it is the simple solution x = A\b
        sol.x(sys.pRCM,1) = sys.A(sys.pRCM,sys.pRCM)\sys.b(sys.pRCM);

        if  exportLinearSol
            if (it==1 && k==1); sol.dx = zeros(size(sys.Al,2),1);end
            bll                 = sys.Al*sol.dx+sys.bl;
            sol.dx(sys.pRCM,1)  = sys.A(sys.pRCM,sys.pRCM)\bll(sys.pRCM);
        end
    end
end