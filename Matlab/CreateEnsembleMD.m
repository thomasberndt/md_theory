function [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V, VBark, HK0, Tc, tau0, N, f, mx, my)
% Creates vectors (V, VBark, HK0, Tc, tau0, N) that describe an ensemble of
% MD particles, based on the same parameters passed into the function. The
% parameters passed into the function can be either vectors or scalars,
% but the output will all be vectors of the right dimensions.
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
% OPTIONAL: 
% f - distribution of grains (scalar or vector) [dimensionless]. Default is 1. 
% mx - normalized magnetization of grains (scalar or vector)
% [dimensionless]. Default is 0.  
% my - additional component for normalized magnetization of grains (scalar or vector)
% [dimensionless]. Default is 0. 
% 
% OUTPUT: 
% The same as above, but all are vectors. 

    if nargin < 9
        my = 0;
    end
    if nargin < 8
        mx = 0;
    end
    if nargin < 7
        f = 1;
    end

    n = max([numel(V), numel(VBark), numel(HK0), numel(Tc), ...
            numel(tau0), numel(N), numel(f), numel(mx), numel(my)]);

    if numel(V)<n
        V = V * ones(1, n); 
    end
    if numel(VBark)<n
        VBark = VBark * ones(1, n); 
    end
    if numel(HK0)<n
        HK0 = HK0 * ones(1, n); 
    end
    if numel(Tc)<n
        Tc = Tc * ones(1, n); 
    end
    if numel(tau0)<n
        tau0 = tau0 * ones(1, n); 
    end
    if numel(N)<n
        N = N * ones(1, n); 
    end
    if numel(f)<n
        f = f * ones(1, n); 
    end
    if numel(mx)<n
        mx = mx * ones(1, n); 
    end
    if numel(my)<n
        my = my * ones(1, n); 
    end
    
    V = V(:)';
    VBark = VBark(:)';
    HK0 = HK0(:)';
    Tc = Tc(:)';
    tau0 = tau0(:)';
    N = N(:)'; 
    f = f(:)';
    mx = mx(:)';
    my = my(:)';
    
end