function mr_out = AcquireTRMMD(V, VBark, HK0, Tc, tau0, N, mr, toT, fromT, t, muH0)
% Simulates TRM acquisition of an ensemble of MD grains given by (V, VBark,
% HK0, Tc, tau0, N). 
% 
% V - volumes of ensemble (vector) [m3]
% VBark - Barkhausen volumes of ensembe (volume swept out by one domainwall
% jump) (vector) [m3]
% HK0 - Microscopic coercivity of pinning at room temperature (vector) [A/m]
% Tc - Curie temperature (vector), used to calculate Ms0 based on
% fitting of data by Dunlop for Titanomagnetite with Ti-content x [K]
% tau0 - attempt time (vector) [s]
% N - shape anisotropy factor (demagnetizing factor) of the domain (vector)
% [dimensionless]
%
% mr - initial state (normalized remanent magnetization) [dimensionless]
% (vector)
%
% toT, fromT - simulates a cooling from fromT to toT [K]
% t - total time of the cooling (linear cooling) [s]
% muH0 - external field (scalar) in [T]
% 
% OUTPUT:
% mr_out - normalized remanent magnetization after TRM acquisition of the
% ensemble (vector) [dimensionless]

    TT = linspace(fromT, toT, 10*min(50, 1+ceil(abs(fromT-toT)/5)));
    is_cooling = (toT < fromT); 
    mr_out = mr;
    TB = BlockingTemperatureMD(V, VBark, HK0, Tc, tau0, N, t./length(TT), mr_out, muH0); 
    
    for n = 1:length(TT)
        % TRM is simulated by a forward model, decreasing the
        % temperature in small steps and calculating the VRM
        % acquisition at each step
        mr_new = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mr_out, TT(n), t./length(TT), muH0);  
        
        if is_cooling && n > 1
            % During the temperature step when the grain is not SP anymore,
            % we need to make sure that its blocked magnetization state
            % corresponds exactly to the the temperature where blocking
            % occured. This calculation must be much more precise than the
            % temperature steps in TT would allow, because the equilibrium
            % magnetization is strongly dependent on temperature, and some
            % grains move from totally unblocked to totally blocked in a
            % fraction of a degree.
            just_blocked = logical(TT(n)<TB & TT(n-1)>TB); 
            mr_eq = EquilibriumMagnetizationMD(V, VBark, Tc, N, TB, muH0); 
            mr_new(just_blocked) = mr_eq(just_blocked); 
        end
        mr_out = mr_new; 
    end
end