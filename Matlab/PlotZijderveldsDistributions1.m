
mu0 = pi*4e-7;
yr = 3600*24*365.25; % Year [s]
color_vrm = [0 0.3 0.7];
color_trm = [0.7 0.3 0];

clf
p = panel();
p.pack(4, 3);
p.marginleft = 24;
p.marginright = 12;
p.margintop = 12; 
p.marginbottom = 18; 
p.de.margin = 12;
p.de.marginleft = 8;
p.de.marginright = 8;
set(gcf, 'Color', 'w');

HK0 = logspace(log10(0.1), log10(10), 20)*1e-3/mu0;
vol = logspace(log10(1e-6^3), log10(100e-6^3), 10);

[HK0, vol] = meshgrid(HK0, vol); 
HK0 = HK0(:);
vol = vol(:); 
f = 1./vol./HK0; 


VBark_set = [ 50e-9^3 70e-9^3 100e-9^3]; 

for n = 1:length(VBark_set)
    p(1, n).select();
    [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(vol, VBark_set(n), HK0, 580 + 273, 1e-10, 0.127, f); 
    [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceVRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 273, 680e3*yr);
    Zijderveld(M, T-273, color_vrm); 
    hold on
    title(sprintf('V_{Bark} = (%g nm)^3, V = (1-100 \\mu{}m)^3', (VBark_set(n)^(1/3))*1e9)); 
    p(3, n).select();
    Zijderveld(Maf, muHaf*1e3, color_vrm);
    hold on
    title(sprintf('V_{Bark} = (%g nm)^3, V = (1-100 \\mu{}m)^3', (VBark_set(n)^(1/3))*1e9)); 

    [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(vol, VBark_set(n), HK0, 580 + 273, 1e-10, 0.127, f); 
    [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceTRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 100+273, 100e3*yr);
    p(1, n).select();
    Zijderveld(M, T-273, color_trm); 
    p(3, n).select();
    Zijderveld(Maf, muHaf*1e3, color_trm);
    drawnow
end



HK0 = logspace(log10(0.1), log10(10), 20)*1e-3/mu0;
volBark = logspace(log10(50e-9^3), log10(100e-9^3), 10);

[HK0, volBark] = meshgrid(HK0, volBark); 
HK0 = HK0(:);
volBark = volBark(:); 
f = 1./volBark./HK0; 


V_set = [1e-6^3 10e-6^3 100e-6^3]; 

for n = 1:length(V_set)
    p(2, n).select();
    [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V_set(n), volBark, HK0, 580 + 273, 1e-10, 0.127, f); 
    [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceVRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 273, 680e3*yr);
    Zijderveld(M, T-273, color_vrm); 
    hold on
    title(sprintf('V_{Bark} = (50-100 nm)^3, V = (%g \\mu{}m)^3', (V_set(n)^(1/3))*1e6)); 
    p(4, n).select();
    Zijderveld(Maf, muHaf*1e3, color_vrm);
    hold on
    title(sprintf('V_{Bark} = (50-100 nm)^3, V = (%g \\mu{}m)^3', (V_set(n)^(1/3))*1e6)); 

    [V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V_set(n), volBark, HK0, 580 + 273, 1e-10, 0.127, f); 
    [M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceTRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 100+273, 100e3*yr);
    p(2, n).select();
    Zijderveld(M, T-273, color_trm); 
    p(4, n).select();
    Zijderveld(Maf, muHaf*1e3, color_trm);
    drawnow
end







try
    export_fig('../output/png/ZijderveldsDistributions1.png', '-m4'); 
    export_fig('../output/pdf/ZijderveldsDistributions1.pdf'); 
catch
    print(gcf, '-dpng', '../output/png/ZijderveldsDistributions1.png'); 
    print(gcf, '-dpdf', '../output/pdf/ZijderveldsDistributions1.pdf'); 
end

