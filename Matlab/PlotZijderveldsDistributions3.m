
mu0 = pi*4e-7;
yr = 3600*24*365.25; % Year [s]
color_vrm = [0 0.3 0.7];
color_trm = [0.7 0.3 0];

clf
p = panel();
p.pack(1, 2);
p.marginleft = 24;
p.marginright = 12;
p.margintop = 12; 
p.marginbottom = 18; 
p.de.margin = 12;
p.de.marginleft = 8;
p.de.marginright = 8;
set(gcf, 'Color', 'w');

HK0 = logspace(log10(0.1), log10(10), 10)*1e-3/mu0;
VBark = logspace(log10(50e-9^3), log10(100e-9^3), 10);
V = logspace(log10(1e-6^3), log10(100e-6^3), 10);

[HK0, VBark, V] = meshgrid(HK0, VBark, V); 
HK0 = HK0(:);
VBark = VBark(:); 
V = V(:); 

f = 1./VBark./HK0./V; 

[V, VBark, HK0, Tc, tau0, N, f, mx, my] = CreateEnsembleMD(V, VBark, HK0, 580 + 273, 1e-10, 0.127, f); 
[M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceVRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 273, 680e3*yr);
p(1, 1).select();
h_t_vrm = Zijderveld(M, T-273, color_vrm); 
ZijderveldLabel(M, T-273, T-273, color_vrm, 'topleft'); 
hold on
title('a) Thermal: V_{Bark} = (50-100nm)^3, V = (1-100\mu{}m)^3'); 
p(1, 2).select();
h_af_vrm = Zijderveld(Maf, muHaf*1e3, color_vrm);
ZijderveldLabel(Maf, muHaf*1e3, muHaf*1e3, color_vrm, 'topleft'); 
title('b) AF: V_{Bark} = (50-100nm)^3, V = (1-100\mu{}m)^3'); 
hold on

[M, T, Maf, muHaf] = PrimaryAndSecondaryRemanenceTRM(V, VBark, HK0, Tc, tau0, N, f, mx, my, 200+273, 100e3*yr);
p(1, 1).select();
h_t_trm = Zijderveld(M, T-273, color_trm); 
ZijderveldLabel(M, T-273, T-273, color_trm, 'bottomright'); 
legend([h_t_vrm, h_t_trm], 'VRM', 'TRM');
p(1, 2).select();
h_af_trm = Zijderveld(Maf, muHaf*1e3, color_trm);
ZijderveldLabel(Maf, muHaf*1e3, muHaf*1e3, color_trm, 'bottomright'); 
legend([h_af_vrm, h_af_trm], 'VRM', 'TRM');
drawnow









try
    export_fig('../output/png/ZijderveldsDistributions3.png', '-m4'); 
    export_fig('../output/pdf/ZijderveldsDistributions3.pdf'); 
catch
    print(gcf, '-dpng', '../output/png/ZijderveldsDistributions3.png'); 
    print(gcf, '-dpdf', '../output/pdf/ZijderveldsDistributions3.pdf'); 
end

