function [Tcrit, tcrit] = CalculateApproximationLimits(V, VBark, HK0, Tc, tau0, N)
% Calculates the critical values where the assumptions made for some of the 
% approximations in the theory do not hold anymore, for an ensemble of 
% MD grains given by (V, VBark, HK0, Tc, tau0, N). 
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
% OUTPUT:
% Tcrit - There is a critical temperature Tcrit [K], 
% above which the approximations are not
% valid any more (because of the deminishing demagnetizing field). Hence, 
% the theory is only valid for T<Tcrit. (vector).
% tcrit - The relaxation time tcrit [s] corresponding to this temperature may be long,
% i.e. corresponding to a blocked grain. This would mean that such a grain
% would be permanently blocked. Hence the theory is not valid at all for
% such a grain: tcrit>a few second => the whole grain is nonsense. 

    mu0 = pi*4e-7; 
    kB = 1.38e-23; 
    Ms0 = CalculateMs0(Tc);
    
    Tcrit = 1./(kB./mu0./HK0./VBark./Ms0 + 1./Tc); 
    tcrit = tau0.*V.*HK0./4./N./Ms0./VBark.*exp(1); 
end