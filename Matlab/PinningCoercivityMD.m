function HK0 = PinningCoercivityMD(V, VBark, Tc, tau0, N, T, t, mr, muH0)
% Calculates what pinning coercivity HK0 would be necessary to obtain the
% given relaxation time and blocking temperatures for the MD
% grain with the given parameters.
% 
% V - volumes of ensemble (vector) [m3]
% VBark - Barkhausen volumes of ensembe (volume swept out by one domainwall
% jump) (vector) [m3]
% Tc - Curie temperature (vector), used to calculate Ms0 based on
% fitting of data by Dunlop for Titanomagnetite with Ti-content x [K]
% tau0 - attempt time (vector) [s]
% N - shape anisotropy factor (demagnetizing factor) of the domain (vector)
% [dimensionless]
%
% T - blocking temperature (scalar or vector) [K]
% t - relaxation time (scalar or vector) [s]
% 
% OPTIONAL PARAMETERS:
% mr - normalized remanent state of the grain (scalar) [dimensionless]
% muH0 - external field (scalar) in [T]
% In non-zero field, there is a dependence on the remanence state. Hence
% either both mr and muH0 need to be given, or non of them.
% 
% OUTPUT:
% HK0 - Microscopic coercivity of pinning at room temperature (vector) [A/m]
    HK00 = logspace(log10(0.001), log(1000), 1000); 
    
    if numel(T) == 1
        T = T * ones(size(t));
    end
    if numel(t) == 1
        t = t * ones(size(T));
    end
        
    HK0 = zeros(size(T)); 
    
    for n = 1:length(T)
        if nargin == 9
            tt = RelaxationTimeMD(V, VBark, HK00, Tc, tau0, N, T(n), mr, muH0);
        else
            tt = RelaxationTimeMD(V, VBark, HK00, Tc, tau0, N, T(n));
        end
        
        idx = find(tt>=t(n),1,'first'); 
        if sum(idx) > 0
            HK0(n) = HK00(idx); 
        end
    end
    
end