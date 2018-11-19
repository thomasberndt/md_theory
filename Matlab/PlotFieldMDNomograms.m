
T = 300;        % Temperature [K]
kB = 1.38e-23;  % Boltzman constant [J/K]
mu0 = pi*4e-7;  % vacuum permeability [T*m/A]
yr = 3600*24*365.25; % Year [s]
Tsteps = (50:50:550) + 273;
t_min = 60;

color_small = repmat([0 0 1], length(Tsteps), 1);
color_large = repmat([0 0.5 0], length(Tsteps), 1);
color_small_bark = repmat([0.7 0.7 0], length(Tsteps), 1);
color_large_bark = repmat([0.5 0 0], length(Tsteps), 1);

clf

p = panel();
p.pack(3, 2);
p.marginleft = 24;
p.marginright = 12;
p.margintop = 12; 
p.marginbottom = 18; 
p.de.margin = 12;
p.de.marginleft = 8;
p.de.marginright = 8;
set(gcf, 'Color', 'w');




VBark_set = [ 50e-9^3 70e-9^3 100e-9^3]; 

for n = 1:length(VBark_set)
    p(n, 1).select();
    hSD = NomogramSD(580 + 273, 1e-10, [0 0 0], ':'); 
    hold on

    [V, VBark, ~, Tc, tau0, N] = CreateEnsembleMD(1e-6^3, VBark_set(n), zeros(size(Tsteps)), 580 + 273, 1e-10, 0.127); 
    HK0 = PinningCoercivityMD(V(1), VBark(1), Tc(1), tau0(1), N(1), Tsteps, t_min);
    h1 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_small);
    h3 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_small, 0, 100e-6);
    
    NomogramLabelsMD(V, VBark, HK0, Tc, tau0, N, num2cell(num2str(round(HK0'*mu0*1000,2)),2), 0.6, color_small); 

    [V, VBark, ~, Tc, tau0, N] = CreateEnsembleMD(100e-6^3, VBark_set(n), zeros(size(Tsteps)), 580 + 273, 1e-10, 0.127); 
    HK0 = PinningCoercivityMD(V(1), VBark(1), Tc(1), tau0(1), N(1), Tsteps, t_min);
    h2 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_large);
    h4 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_large, 0, 100e-6);
    NomogramLabelsMD(V, VBark, HK0, Tc, tau0, N, num2cell(num2str(round(HK0'*mu0*1000,2)),2), 0.4, color_large); 
    title(sprintf('V_{Bark} = %g nm', VBark_set(n)^(1/3)*1e9)); 
    xlabel('');
    legend([hSD(1), h1(1), h2(1)], {'SD', 'V = 1 \mu{}m', 'V = 100 \mu{}m'}); 
    hold off
    drawnow
end




V_set = [1e-6^3 10e-6^3 100e-6^3]; 

for n = 1:length(V_set)
    p(n, 2).select();
    hSD = NomogramSD(580 + 273, 1e-10, [0 0 0], ':'); 
    hold on

    [V, VBark, ~, Tc, tau0, N] = CreateEnsembleMD(V_set(n), 50e-9.^3, zeros(size(Tsteps)), 580 + 273, 1e-10, 0.127); 
    HK0 = PinningCoercivityMD(V(1), VBark(1), Tc(1), tau0(1), N(1), Tsteps, t_min);
    h1 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_small_bark);
    h3 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_small_bark, 0, 100e-6);
    NomogramLabelsMD(V, VBark, HK0, Tc, tau0, N, num2cell(num2str(round(HK0'*mu0*1000,2)),2), 0.6, color_small_bark); 

    [V, VBark, ~, Tc, tau0, N] = CreateEnsembleMD(V_set(n), 100e-9.^3, zeros(size(Tsteps)), 580 + 273, 1e-10, 0.127); 
    HK0 = PinningCoercivityMD(V(1), VBark(1), Tc(1), tau0(1), N(1), Tsteps, t_min);
    h2 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_large_bark);
    h4 = NomogramMD(V, VBark, HK0, Tc, tau0, N, color_large_bark, 0, 100e-6);
    NomogramLabelsMD(V, VBark, HK0, Tc, tau0, N, num2cell(num2str(round(HK0'*mu0*1000,2)),2), 0.4, color_large_bark); 
    title(sprintf('V = %g \\mu{}m', V_set(n)^(1/3)*1e6)); 
    xlabel('');
    ylabel('');
    legend([hSD(1), h1(1), h2(1)], {'SD', 'V_{Bark} = 50 nm', 'V_{Bark} = 100 nm'}); 
    hold off
    drawnow
end






try
    export_fig('../output/png/Nomograms.png', '-m4'); 
    export_fig('../output/pdf/Nomograms.pdf'); 
catch
    print(gcf, '-dpng', '../output/png/Nomograms.png'); 
    print(gcf, '-dpdf', '../output/pdf/Nomograms.pdf'); 
end




