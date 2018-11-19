function T = BlockingTemperatureMD(V, VBark, HK0, Tc, tau0, N, t, mr, muH0)
% Calculates the blocking temperature of an ensemble of MD grains given by (V, VBark,
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
% t - relaxation time [s]
% 
% OPTIONAL PARAMETERS:
% mr - normalized remanent state of the ensemble (vector) [dimensionless]
% muH0 - external field (scalar) in [T]
% In non-zero field, there is a dependence on the remanence state. Hence
% either both mr and muH0 need to be given, or non of them.
% 
% OUTPUT:
% T - blocking temperatures of the ensemble (vector) [K]    
    TT = linspace(0, max(Tc), 10*max(Tc)); 
    T = zeros(size(V)); 
        
    for n = 1:length(V)
        if nargin == 9
            tt = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), TT, mr(n), muH0);
        else
            tt = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), TT);
        end
        
        idx = find(tt<=t,1,'first'); 
        if sum(idx) > 0
            T(n) = TT(idx); 
        end
    end
    
end