function turbInputSet = controlSet_greedy(Wp)
    turbInputSet = struct();
    
    nTurbs = length(Wp.turbine.Crx);
    turbInputSet.t = 0:Wp.sim.h:1e4; % 10,000 seconds of simulation
    turbInputSet.phi = zeros(nTurbs,length(turbInputSet.t));
    turbInputSet.CT_prime = 2.0*ones(nTurbs,length(turbInputSet.t));
    turbInputSet.interpMethod = 'nearest'; % Nearest since no interpolation necessary
end