function [Mthermal, Tsteps, Maf, muHaf] = PrimaryAndSecondaryRemanenceVRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, T, t, muH0, Tsteps, muHaf)
% Simulates a primary TRM, followed by a perpendicular secondary VRM for
% the ensemble of MD particles. Then returns the stepwise thermal
% demagnetization and the AF demagnetization for it. 
%
% V - volumes of ensemble (scalar or vector) [m3]
% VBark - Barkhausen volumes of ensembe (volume swept out by one domainwall
% jump) (scalar or vector) [m3]
% HK0 - Microscopic coercivity of pinning at room temperature (vector) [A/m]
% Tc - Curie temperature (scalar or vector), used to calculate Ms0 based on
% fitting of data by Dunlop for Titanomagnetite with Ti-content x [K]
% tau0 - attempt time (scalar or vector) [s]
% N - shape anisotropy factor (demagnetizing factor) of the domain (scalar or vector)
% [dimensionless]
% 
% f - grain distribution (vector) [dimensionless] 
% mx, my - normalized remanences of the ensemble (vectors).
% 
% T - temperature for secondary remanence acquisition (scalar) [K]
% t - time it takes to cool from T to room temperature during secondary
% remanence acquisition (scalar) [s]
%
% OPTIONAL:
% muH0 - external field (scalar) in [T]. Default is 30 uT.
% Tsteps - temperature steps of the stepwise thermal demagnetization
% (vector) [K]. Default is 50 K steps. 
% muHaf - AF peak fields (vector) [T] of the AF demagnetization
% (vector) [K]. Default is 0-10 mT in 1 mT steps. 
% 
% OUTPUT: 
% Mthermal - magnetization vectors of stepwise thermal demagnetization
% experiment [Am2]
% Tsteps - temperature steps of the stepwise thermal demagnetization
% (vector) [K]
% Maf- magnetization vectors of AF demagnetization experiment [Am2]
% muHaf - AF peak fields (vector) [T] of the AF demagnetization. 

    if nargin < 14
        muHaf = [0 0.5 1 1.5 2:1:10] * 1e-3;
    end
    if nargin < 13
        Tsteps = (0:25:580) + 273; 
    end
    if nargin < 12
        muH0 = 30e-6; 
    end

    T0 = 273;        % Temperature [K]
    heating_rate = 1; 
    t_hold = 300; 

    primary_cooling_time = 3600; 

    mx = AcquireTRMMD(V, VBark, HK0, Tc, tau0, N, mx, T0, max(Tc), primary_cooling_time, muH0); 
    my = AcquireTRMMD(V, VBark, HK0, Tc, tau0, N, my, T0, max(Tc), primary_cooling_time, 0); 

    mx = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mx, T, t, 0); 
    my = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, my, T, t, muH0); 

    Mx = StepwiseDemagnetization(V, VBark, HK0, Tc, tau0, N, mx, f, T0, Tsteps, heating_rate, t_hold);
    My = StepwiseDemagnetization(V, VBark, HK0, Tc, tau0, N, my, f, T0, Tsteps, heating_rate, t_hold);
    
    Mx_af = AfDemagnetization(V, VBark, HK0, Tc, tau0, N, mx, f, T0, muHaf); 
    My_af = AfDemagnetization(V, VBark, HK0, Tc, tau0, N, my, f, T0, muHaf); 
    
    Mthermal = [Mx', My']; 
    Maf = [Mx_af', My_af']; 
end