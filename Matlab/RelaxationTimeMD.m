function t = RelaxationTimeMD(V, VBark, HK0, Tc, tau0, N, T, mr, muH0)
% Calculates the relaxation time of an ensemble of MD grains given by (V, VBark,
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
% T - temperature [K]
% 
% OPTIONAL PARAMETERS:
% mr - normalized remanent state of the ensemble (vector) [dimensionless]
% muH0 - external field (scalar) in [T]
% In non-zero field, there is a dependence on the remanence state. Hence
% either both mr and muH0 need to be given, or non of them.
% 
% OUTPUT:
% t - relaxation times of the ensemble (vector) [s]

    mu0 = pi*4e-7; 
    kB = 1.38e-23; 
    Ms0 = CalculateMs0(Tc); 
    [Ms, beta] = CalculateMsT(T, Tc);

    if nargin == 9 && muH0~=0
        % In field

        H0 = muH0 ./ mu0;
        d = mu0.*N.*Ms.^2.*VBark.^2 ./ (kB.*T.*V); 
        c = (1+2.*d)./(1+4.*d); 

        t_relax1 = tau0 .* kB .* T .* V ./  ...
            (4 .* mu0 .* N .* Ms0.^2 .* beta.^2 .* VBark.^2) .* ...
            exp(mu0 .* VBark .* Ms0 .* (HK0 .* beta.^2 + ...
            N.*VBark.*Ms./(2.*V)) ./ (kB .* T) );
        fieldrelax = sech(...
            -mu0 .* VBark .* Ms .* (N.*mr - H0) ./ ...
            (kB .* T));
        t_relax2 = t_relax1 .* fieldrelax;
        if size(beta) == size(t_relax2)
            t_relax2(beta==0) = 0;
        else
            if beta == 0
                t_relax2 = zeros(size(t_relax2)); 
            end
        end
        t_relax3 = t_relax2 ./ (1+4.*d); 
        t_relax3(isnan(t_relax3)) = Inf;
        t = t_relax3; 
    else
        % Zero-field

        t = tau0 .* kB .* T .* V ./  ...
                (4 .* mu0 .* N .* Ms.^2 .* VBark.^2) .* ...
                exp(mu0 .* VBark .* Ms0 .* HK0 .* beta.^2 ./ (kB .* T));
    end
    
    % Above some critical temperature Tcrit, some of the approximations
    % made in the MD theory do not hold anymore. The relaxation times above
    % this critical value are set to the smallest relaxation time where the
    % approximations are still valid (which should correspond to
    % superparamagnetic grains anyways - if not, the theory should not be
    % used for this grain). 
    [Tcrit, tcrit] = CalculateApproximationLimits(V, VBark, HK0, Tc, tau0, N);
    if numel(tcrit) > 1
        t(T>Tcrit) = tcrit(T>Tcrit); 
    else
        t(T>Tcrit) = tcrit; 
    end
end


