function h = NomogramSD(Tc, tau0, color, linestyle) 
% Plots a nomogram an ensemble of SD grains with the given Tc and tau0.
%
% Tc - Curie temperature (scalar), used to calculate Ms0 based on
% fitting of data by Dunlop for Titanomagnetite with Ti-content x [K]
% tau0 - attempt time (scalar) [s]
% 
% OPTIONAL:
% color - 3-element vector giving the color for the lines
% 
% OUTPUT:
% h - vector containing the handles of all the lines.
    
    ho = ishold(); 

    yr = 365.25*24*3600;
    T = linspace(0,Tc);
    T_lines = [50:50:550]+273; 
    t_min = 60;
    
    h = zeros(size(T_lines)); 
    
    for n = 1:length(T_lines)
        c = log(t_min/tau0) / (1/T_lines(n) - 1/Tc); 
        lnt = log(tau0) + c * (1./T-1/Tc);
        if nargin >= 3 
            mycolor = color;
        else
            mycolor = [0 0 0];
        end
        if nargin < 4 
            linestyle = '-'; 
        end
        h(n) = semilogy(T-273, exp(lnt), linestyle,'Color', mycolor, 'LineWidth', 1);
        hold on
    end
    
    
    set(gca,'YTick',[1,60,3600,24*3600, 30*24*3600, ...
            yr,10*yr, 100*yr, 1000*yr, 1e4*yr, 1e5*yr, 1e6*yr, ...
            1e7*yr, 1e8*yr 1e9*yr, 4.6e9*yr]);
    set(gca,'YTickLabel',{'1 s', '1 min', '1 hour', '1 day', '1 month', ...
            '1 yr', '10 yr', '100 yr', ...
            '1 ka', '10 ka', '100 ka', ...
            '1 Ma', '10 Ma', '100 Ma', '1 Ba', '4.6 Ba'});

    axis([0,580,60,4.6e9*yr]);
    set(gca, 'YScale', 'log'); 
    
    xlabel('Blocking Temperature T_B [C]'); 
    ylabel('Relaxation time t [s]'); 
    
    if ho
        hold on
    else
        hold off
    end
end