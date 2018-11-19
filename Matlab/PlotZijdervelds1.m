
mu0 = pi*4e-7;
yr = 3600*24*365.25; % Year [s]
color_vrm = [0 0.3 0.7];
color_trm = [0.7 0.3 0];

f1 = figure(1); 
clf
p1 = panel();
p1.pack(3, 3);
p1.marginleft = 24;
p1.marginright = 12;
p1.margintop = 12; 
p1.marginbottom = 18; 
p1.de.margin = 12;
p1.de.marginleft = 8;
p1.de.marginright = 8;
set(gcf, 'Color', 'w');

f2 = figure(2); 
clf
p2 = panel();
p2.pack(3, 3);
p2.marginleft = 24;
p2.marginright = 12;
p2.margintop = 12; 
p2.marginbottom = 18; 
p2.de.margin = 12;
p2.de.marginleft = 8;
p2.de.marginright = 8;
set(gcf, 'Color', 'w')


VBark_set = [ 50e-9^3 70e-9^3 100e-9^3]; 
V_set = [1e-6^3 10e-6^3 100e-6^3]; 
HK0 = logspace(log10(0.1), log10(6))*1e-3/mu0;

letters = 'abcdefghi'; 

for n = 1:length(VBark_set)
    for k = 1:length(V_set)
        [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V_set(k), VBark_set(n), HK0, 580 + 273, 1e-10, 0.127); 
        [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceVRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 273, 680e3*yr);
        p1(k, n).select();
        h_t_vrm = Zijderveld(M, T-273, color_vrm); 
        ZijderveldLabel(M, T-273, T-273, color_vrm, 'topleft');
        hold on
        title(sprintf('%s) Thermal: V_{Bark} = (%g nm)^3, V = (%g \\mu{}m)^3', ...
            letters(3*k+n-3), (VBark_set(n)^(1/3))*1e9, (V_set(k)^(1/3))*1e6)); 
        p2(k, n).select();
        h_af_vrm = Zijderveld(Maf, muHaf*1e3, color_vrm);
        ZijderveldLabel(Maf, muHaf*1e3, muHaf*1e3, color_vrm, 'topleft');
        hold on
        title(sprintf('%s) AF: V_{Bark} = (%g nm)^3, V = (%g \\mu{}m)^3', ...
            letters(3*k+n-3), (VBark_set(n)^(1/3))*1e9, (V_set(k)^(1/3))*1e6)); 
           
        [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V_set(k), VBark_set(n), HK0, 580 + 273, 1e-10, 0.127); 
        [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceTRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 100+273, 100e3*yr);
        p1(k, n).select();
        h_t_trm = Zijderveld(M, T-273, color_trm); 
        ZijderveldLabel(M, T-273, T-273, color_trm);
        legend([h_t_vrm, h_t_trm], 'VRM', 'TRM');
        p2(k, n).select();
        h_af_trm = Zijderveld(Maf, muHaf*1e3, color_trm);
        ZijderveldLabel(Maf, muHaf*1e3, muHaf*1e3, color_trm);
        legend([h_af_vrm, h_af_trm], 'VRM', 'TRM');
        drawnow
    end
end


figure(f1);
try
    export_fig('../output/png/ZijderveldsThermal.png', '-m4'); 
    export_fig('../output/pdf/ZijderveldsThermal.pdf'); 
catch
    print(gcf, '-dpng', '../output/png/ZijderveldsThermal.png'); 
    print(gcf, '-dpdf', '../output/pdf/ZijderveldsThermal.pdf'); 
end

figure(f2);
try
    export_fig('../output/png/ZijderveldsAF.png', '-m4'); 
    export_fig('../output/pdf/ZijderveldsAF.pdf'); 
catch
    print(gcf, '-dpng', '../output/png/ZijderveldsAF.png'); 
    print(gcf, '-dpdf', '../output/pdf/ZijderveldsAF.pdf'); 
end


