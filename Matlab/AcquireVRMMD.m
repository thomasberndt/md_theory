function mr_out = AcquireVRMMD(V, VBark, HK0, Tc, tau0, N, mr, T, t, muH0)
% Simulates VRM acquisition of an ensemble of MD grains given by (V, VBark,
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
% T - temperature during VRM acquisition [K]
% t - total time of VRM acquisition [s]
% muH0 - external field (scalar) in [T]
% 
% OUTPUT:
% mr_out - normalized remanent magnetization after VRM acquisition of the
% ensemble (vector) [dimensionless]
    mu0 = pi*4e-7; 
    kB = 1.38e-23;
    [ms] = CalculateMsT(T, Tc);
    H0 = muH0 ./ mu0;
    
    d = mu0.*N.*ms.^2.*VBark.^2 ./ (kB.*T.*V); 
    c = (1+2.*d)./(1+4.*d); 
    t_relax = RelaxationTimeMD(V, VBark, HK0, Tc, tau0, N, T, mr, muH0);

    mr_eq = ones(size(mr)).* c .* H0 ./ N ./ ms; 
    mr_eq(ms==0) = 0;
    saturated = logical(abs(mr_eq) >= 1); 
    mr_eq(saturated) = sign(mr_eq(saturated));
    
    relax = exp(-t./t_relax); 
    relax(t_relax==0) = 0;
    mr_out = mr .* relax + mr_eq .* (1-relax);
    saturated = logical(abs(mr_out) >= ms); 
    mr_out(saturated) = ms(saturated) .* sign(mr_out(saturated));
    mr_out = real(mr_out); 
end