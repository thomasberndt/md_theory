function Mr = AfDemagnetization(V, VBark, HK0, Tc, tau0, N, mr, f, T, muH)
% Simulates AF demagnetization of an ensemble of MD grains 
% given by (V, VBark, HK0, Tc, tau0, N). Thermal relaxation is neglected
% for AF demagnetization. 
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
% T - the temperature at which the remanence is measured (i.e. room 
% temperature) (scalar) [K] 
% muH - the AF peak fields (vector) [T]
% 
% OUTPUT:
% Mr - total remanent magnetization vectors after each demagnetization step
% (vector) [Am2]

    mu0 = pi*4e-7; 
    t_measure = 100; % time it takes to measure the sample

    H = muH/mu0; 
    Mr = zeros(size(muH));  
    
    % Viscous decay during the time it takes to put sample into the
    % instrument
    mr = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mr, T, t_measure, 0);  
    
    for n = 1:length(muH)
        remanence = (HK0 > H(n)); 
        Mr(n) = MeasureNRM(mr.*remanence, Tc, f, V);  
    end
end