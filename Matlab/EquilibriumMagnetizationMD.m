function mr_eq = EquilibriumMagnetizationMD(V, VBark, Tc, N, T, muH0)
% Calculates the Equilibrium magnetization of an ensemble of MD grains 
% given by (V, VBark, Tc, N) (the equilibrium magnetization does not depend 
% on HK0 and tau0). 
% Note that for small fields, the equilibrium magnetization is
% approximately H0/N. 
% 
% V - volumes of ensemble (vector) [m3]
% VBark - Barkhausen volumes of ensembe (volume swept out by one domainwall
% jump) (vector) [m3]
% Tc - Curie temperature (vector), used to calculate Ms0 based on
% fitting of data by Dunlop for Titanomagnetite with Ti-content x [K]
% N - shape anisotropy factor (demagnetizing factor) of the domains (vector)
% [dimensionless]
%
% T - temperature [K]
% muH0 - external field [T]
%
% OUTPUT: 
% mr_eq - the normalized equilibrium magnetization (vector) [dimensionless]
    mu0 = pi*4e-7; 
    kB = 1.38e-23;
    [ms, ~] = CalculateMsT(T, Tc);
    H0 = muH0 ./ mu0;
    
    d = mu0.*N.*ms.^2.*VBark.^2 ./ (kB.*T.*V); 
    c = (1+2.*d)./(1+4.*d); 
    
    mr_eq = ones(size(V)).* c .* H0 ./ N ./ ms; 
end