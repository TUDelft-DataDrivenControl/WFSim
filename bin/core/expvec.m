function [ Ld ] = expvec( L, cellSize, Ct, R_con, N_con, dx_con_min )
%EXPVEC Calculate quasi-unidistant vector with exponential concentrations
%   This function determines a discrete number of points in a domain with
%   length 'L', with equidistant spacing 'cellSize'. Points are
%   concentrated around turbines Ct(1)...Ct(n), with concentration radius
%   'R_con', and within this radius a total number of points 'N_con'. The
%   minimum distance between two points is 'dx_con_min' in this radius.
%
%   Note: this function returns a linear grid when Ct = []
%
%   Date: October 18th, 2016
%   Author: Bart
%

if nargin <= 0
    L           = 70;      % Domain length
    cellSize    = 5;       % Approximate cell size
    Ct          = [20,40]; % Turbine locations

    R_con       = 3;       % Concentration radius x-direction
    N_con       = 11;      % Number of concentrated points around in your radius (sum of left and right)
    dx_con_min  = 0.03;    % Minimal spacing around turbine
elseif dx_con_min*(N_con+2) > R_con*2
    error('dx_con_min is too large to create a dense spacing around the turbines.')
elseif sum(diff(Ct) <= 2*R_con) > 0
    error('Please put turbines further apart, at identical locations, or make the concentration radius smaller.');
end;

Ld = 0;
for i = 1:length(Ct)
    % Equidistant spacing from current until the next concentration region
    Ld = [Ld(1:end-1) linspace(Ld(end),Ct(i)-R_con,round((Ct(i)-R_con-Ld(end))/cellSize)+1)];

    % Setup a nonlin function to determine spacing in concentrated region
    syms y(a,b,c,x)
    centerPiece = (N_con-1)/2+1;
    
    y  = symfun(a+b*x+c*x^2,x); %y = symfun(a+b*x+c*x^2+d*x^1.05,x); %y = symfun(a+b*x^c,x);
    yd = diff(y,x);
    p  = solve([y(0) == Ct(i)-R_con, ...
               y(centerPiece)  == Ct(i), ... 
               ... %yd(0) == (ldx(2)-ldx(1)), ...
               yd(centerPiece) == dx_con_min],[a,b,c]);

    y = subs(y, a, double(p.a));
    y = subs(y, b, double(p.b));
    y = subs(y, c, double(p.c));
    %y = subs(y, d, double(p.d));

    % Process datapoints in concentrated region
    ldx_l           = double(y([1:centerPiece]));
    ldx_r           = 2*Ct(i)-ldx_l; 
    ldx_r(end:-1:1) = ldx_r; % swap order
    
    % Merge datapoints into total vector
    if ( round(N_con/2) == (N_con/2)) 
        Ld = [Ld ldx_l ldx_r Ct(i)+R_con];
    else
        Ld = [Ld ldx_l ldx_r(2:end) Ct(i)+R_con];
    end;
end;
% Finally, add equidistant points until end of domain
Ld = [Ld(1:end-1) linspace(Ld(end),L,round((L-Ld(end))/cellSize)+1)];
end