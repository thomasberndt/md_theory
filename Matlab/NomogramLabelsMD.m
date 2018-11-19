function h = NomogramLabelsMD(V, VBark, HK0, Tc, tau0, N, labels, position, color) 
% Adds labels to a nomogram for the ensemble of MD grains given by 
% (V, VBark, HK0, Tc, tau0, N). 
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
% labels - array of strings to label the lines
% 
% OPTIONAL:
% position - value between 0 and 1 indicating where the labels are placed:
% 0 corresponds to the top left end of the lines, 1 to the bottom right end
% of the lines. Default is 0.5.
% color - Tx3 matrix giving the color for each grain in the ensemble
% 
% OUTPUT:
% h - vector containing the handles of all the lines.
    
    gW = 0.93e-3; % energy of a domain wall
    yr = 365.25*24*3600;
    
    if nargin < 8
        position = 0.5;
    end
    
    ho = ishold(); 
    
    h = zeros(size(V));
    
    for n = 1:length(V)
        xl = xlim;
        T0 = xl(1)+273;
        t0 = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), T0);
        yl = ylim; 
        t1 = yl(1);
        t2 = min(t0, yl(2)); 
        T1 = BlockingTemperatureMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), t1);
        T2 = BlockingTemperatureMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), t2);
        Tp = T1 + (T2-T1)*position; 
        dTp = Tp + 1; 
        tp  = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), Tp);
        dtp = RelaxationTimeMD(V(n), VBark(n), HK0(n), Tc(n), tau0(n), N(n), dTp); 
        nor = (xl(2)-xl(1)) / log(yl(2)/yl(1));
        angle = atan(nor*log(tp/dtp));
        if nargin >= 9
            mycolor = color(n,:); 
        else 
            mycolor = [0 0 0]; 
        end
        h(n) = text(Tp-273, tp, labels{n}, 'Color', mycolor, 'Rotation', -angle/pi*180, 'FontSize', 12);
        hold on
    end
    
    if ho
        hold on
    else
        hold off
    end
end