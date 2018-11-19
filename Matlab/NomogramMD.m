function h = NomogramMD(V, VBark, HK0, Tc, tau0, N, color, mr, muH0) 
% Plots a nomogram for the ensemble of MD grains given by (V, VBark,
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
% OPTIONAL:
% color - Tx3 matrix giving the color for each grain in the ensemble
% mr - normalized remanent state of the ensemble (vector) [dimensionless]
% muH0 - external field (scalar) in [T]
% In non-zero field, there is a dependence on the remanence state. Hence
% either both mr and muH0 need to be given, or non of them.
% 
% OUTPUT:
% h - vector containing the handles of all the lines.
    
    gW = 0.93e-3; % energy of a domain wall
    yr = 365.25*24*3600;
    
    ho = ishold(); 
    
    h = zeros(size(V));
        
    for n = 1:length(V)
        T = linspace(0, Tc(n)-1, 1000);
        
        if nargin < 9
            t = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), T);
        else
            t = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), T, mr, muH0);
        end
        if nargin >= 7
            mycolor = color(n,:); 
        else 
            mycolor = [0 0 0]; 
        end
        h(n) = semilogy(T-273, t,'Color', mycolor, 'LineWidth', 1);
        hold on
    end
    
    % Make graphs clickable
    set(h, 'ButtonDownFcn', {@LineSelected, h})
    
    % Use sensible timescale
    set(gca,'YTick',[1,60,3600,24*3600, 30*24*3600, ...
            yr,10*yr, 100*yr, 1000*yr, 1e4*yr, 1e5*yr, 1e6*yr, ...
            1e7*yr, 1e8*yr 1e9*yr, 4.5e9*yr]);
    set(gca,'YTickLabel',{'1 s', '1 min', '1 hour', '1 day', '1 month', ...
            '1 yr', '10 yr', '100 yr', ...
            '1 ka', '10 ka', '100 ka', ...
            '1 Ma', '10 Ma', '100 Ma', '1 Ba', '4.6 Ba'});

%     grid on
    axis([0,580,60,4.5e9*yr]);
    set(gca, 'YScale', 'log'); 
    
    xlabel('Blocking Temperature T_B [C]'); 
    ylabel('Relaxation time t [s]'); 
    
    if ho
        hold on
    else
        hold off
    end
    
    % When a line is clicked, show the correpsponding grain parameters in
    % the title
    function LineSelected(ThisLine, ~, Lines)
        set(ThisLine, 'LineWidth', 2.5);
        set(Lines(Lines ~= ThisLine), 'LineWidth', 1);
        a = find(Lines == ThisLine); 
        title(sprintf('V: %g um, VBark: %g nm, HK0: %g uT, %g', ...
                round(V(a).^(1/3) / 1e-6), ...
                round(VBark(a).^(1/3) / 1e-9), ...
                round(HK0(a) * pi*4e-7 * 1000000)));
    end
    
    
end