function Mr = StepwiseDemagnetization(V, VBark, HK0, Tc, tau0, N, mr, f, T0, T, heating_rate, t, muH0)
% Simulates stepwise thermal demagnetization of an ensemble of MD grains 
% given by (V, VBark, HK0, Tc, tau0, N). 
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
% f - the grain distribution of the grains in the above ensemble (vector)
% [dimensionless]
%
% T0 - the temperature at which the remanence is measured (i.e. room 
% temperature) (scalar) [K] 
% T - the temperature steps (vector) [K] 
% heating_rate - the heating (and cooling) rate of the furnace (scalar)
% [K/s]
% t - hold total time of the temperature (scalar) [s]
% 
% OPTIONAL: 
% muH0 - applied field (scalar) during heating, hold and cooling in [T].
% Default is zero-field
% 
% OUTPUT:
% Mr - total remanent magnetization vectors after each heating step
% (vector) [Am2]

    if nargin < 13 
        muH0 = 0; 
    end

    t_measure = 100; % time it takes to measure the sample
    
    Mr = zeros(size(T));
    
    for n = 1:length(T) 
        dT = T(n) - T0;       % how many degrees from room T up to target T?
        dt = dT/heating_rate; % time it takes to heat up to T from room T
        
        % Heat up
        mr = AcquireTRMMD(V, VBark, HK0, Tc, tau0, N, mr, T(n), T0, dt, muH0);
        % Hold 
        mr = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mr, T(n), t, muH0);
        % Cool down
        mr = AcquireTRMMD(V, VBark, HK0, Tc, tau0, N, mr, T0, T(n), dt, muH0);
        % Put the sample into zero-field
        mr = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mr, T0, t_measure, 0);
        % Measure the sample
        Mr(n) = MeasureNRM(mr, Tc, f, V); 
    end
end