function lm=ConstructLmu(x_IF,y_IF,WD,xTurbs,yTurbs,D,d_lower,d_upper,lm_slope)
    % Check inputs
    if d_upper <= d_lower
        error('Make sure your upper bound is larger than your lower bound on the Lmu turbulence model.');
    end

    lm = zeros(size(x_IF));    
    
    % Rotation vector
    rotWD = [cos(WD) -sin(WD);sin(WD) cos(WD)];
    
    % Add mixing length to the field for every turbine
    for iT = 1:length(xTurbs)
        xy_WF = [x_IF(:)-xTurbs(iT) y_IF(:)-yTurbs(iT)]*rotWD;
        x_WF  = reshape(xy_WF(:,1),size(x_IF,1),size(x_IF,2));
        y_WF  = reshape(xy_WF(:,2),size(y_IF,1),size(y_IF,2));

        % Determine turbine-added mixing length
        lm = lm + Lmu_2D_WF(x_WF,y_WF,D,d_lower,d_upper,lm_slope);
    end
 
    H     = diskfilter('disk',1);
    lm    = filter2(H,lm);
    
    % % Plot results
    % clf; surf(X,Y,lm);
    % axis equal tight
    % xlabel('x (m)');
    % ylabel('y (m)');
    % title('Lmu plot (m)');
    % drawnow()

    % 2D Lmu profile for single turbine case at (x,y)=(0,0)
    %    this is the classical profile with sharp corners
    function lm=Lmu_2D_WF(x,y,D,d_lower,d_upper,lm_slope)
        lm       = zeros(size(x));
        indx     = ((x > d_lower) & (x < d_upper) & (y <= D/2) & (y > -D/2));
        lm(indx) = (x(indx)-d_lower).*lm_slope;
    end

    % 2D Lmu profile for single turbine case at (x,y)=(0,0)
    %    this is the quadratic shape with rounded corners and expanding width
%     function lm=Lmu_2D_WF(x,y,D,d_lower,d_upper,lm_slope)
%         k_2      = 0.003;
%         lm       = zeros(size(x));
%         indx = ((x > d_lower) & (x < d_upper) & (abs(y) < sqrt(((x-d_lower).*lm_slope)/k_2)));
%         lm(indx) = (x(indx)-d_lower)*lm_slope - (y(indx).^2)*k_2;
%     end
end
