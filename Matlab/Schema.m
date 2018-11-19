
figure(1)
clf

p = panel();
p.pack(3, 1);
p.marginleft = 8;
p.marginright = 6;
p.margintop = 6; 
p.marginbottom = 12; 
p.de.margin = 12;
p.de.marginbottom = 18;
set(gcf, 'Color', 'w');


p(1,1).select();


theta = linspace(0, 2*pi,1000); 
Theta = theta / pi * 180-90;

H = 0.05 * theta; 
E = 0.3*cos(2*theta) + H + 0.4; 

[~,ids] = findpeaks(-E); 
E1 = E(ids(1));
E2 = E(ids(2));
T1 = Theta(ids(1));
T2 = Theta(ids(2));

plot(Theta, E, 'k-');
hold on
plot([T1 T2], [E1 E2]+0.04, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 10); 
plot(Theta, H, 'k--');

xlim([min(Theta) max(Theta)]); 
ylim([-0.1*max(E) max(E)*1.1]);

xlabel('Angle from the easy axis [°]'); 
ylabel('Energy [arbitrary units]'); 

text(T1, E1-0.03, 'n_+', 'HorizontalAlignment', 'Center');
text(T2, E2-0.03, 'n_-', 'HorizontalAlignment', 'Center');
text(90, 0.1, 'External field', 'HorizontalAlignment', 'Center');
text(93, 0.6, 'Energy barrier', 'HorizontalAlignment', 'Center');
text(93, 1, 'Flipping of states due to thermal activations', 'HorizontalAlignment', 'Center');

arrowT = linspace(140, 220, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.3*cos(1+1.69*arrowt) + 0.05*arrowt + 0.5; 
aE = arrowE(1); 
aT = arrowT(1)-90; 

plot(arrowT-90, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT+10 aT], [aE aE+0.02 aE+0.05], 'k'); 

ax = gca;
ax.XTick = -90:45:270;
ax.YTick = [];

title('a) Néel SD theory'); 











p(2,1).select();

Eampl = 0.05; 
H = linspace(-0.1, 0.1,200);
E = Eampl*cos(370*H) + 30*H.^2 + 1*H + 0.2; 

[~,ids] = findpeaks(-E); 
N = 7;
E1 = E(ids(N));
E2 = E(ids(N+1));
H1 = H(ids(N));
H2 = H(ids(N+1));

plot(H*10, E, 'k-');
hold on
plot([H1 H2]*10, [E1 E2]+0.03, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 10); 
plot(H*10, H+0.05, 'k--');

xlim([min(H) max(H)]*10); 
ylim([-0.1*max(E) max(E)*1.1]);

xlabel('Distance from centre of the grain [nm]'); 
ylabel('Energy [arbitrary units]'); 

text(H1*10, E1-0.02, 'n_i', 'HorizontalAlignment', 'Center');
text(H2*10, E2-0.02, 'n_{i+1}', 'HorizontalAlignment', 'Center');
text(0, 0.02, 'External field', 'HorizontalAlignment', 'Center');
text(-0.52, 0.14, 'Energy barrier', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.60, 0.31, 'Pinning site', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.69, 0.20, 'Energy barrier', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.77, 0.36, 'Pinning site', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(0.17, 0.47, {'Flipping of states due ', 'to thermal activations'}, ...
            'HorizontalAlignment', 'Center', ...
            'Rotation', 90);

arrowT = linspace(140, 240, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.05*cos(1+1.69*arrowt) + 0.05*arrowt + 0.1; 
aE = arrowE(1); 
aT = (arrowT(1)-140)/800+0.1; 

plot((arrowT-140)/800+0.1, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT+0.04 aT], [aE aE+0.02 aE+0.03], 'k'); 

ax = gca;
ax.XTick = -1:0.5:1;
ax.YTick = [];

title('b) Néel MD theory'); 













p(3,1).select();

H = linspace(-0.1, 0.1,200);
Eampl = 0.05; 
E = Eampl*cos(370*H) + 30*H.^2 + 1*H + 0.2; 

[~,ids] = findpeaks(-E); 
N = 7;
EE = E(ids(N+(-2:2)));
HH = H(ids(N+(-2:2)));
pp = [1 5 15 5 1]; 

plot(H*10, E, 'k-');
hold on
for n = 1:length(EE)
    rE = 1.5 * Eampl * (rand(pp(n), 1)-0.2); 
    rH = (HH(2)-HH(1))/2 * (rand(pp(n), 1)-0.5); 
    if pp(n) == 1
        rE = 0;
        rH = 0;
    end
    plot((HH(n)+rH)*10, EE(n)+rE+0.03, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end
plot(H*10, H+0.05, 'k--');

xlim([min(H) max(H)]*10); 
ylim([-0.1*max(E) max(E)*1.1]);

xlabel('Distance from centre of the grain [nm]'); 
ylabel('Energy [arbitrary units]'); 

text(H(ids(N))*10, E(ids(N))-0.02, 'n_i', 'HorizontalAlignment', 'Center');
text(H(ids(N+1))*10, E(ids(N+1))-0.02, 'n_{i+1}', 'HorizontalAlignment', 'Center');
text(H(ids(N-1))*10, E(ids(N-1))-0.02, 'n_{i-1}', 'HorizontalAlignment', 'Center');
text(0, 0.02, 'External field', 'HorizontalAlignment', 'Center');
text(-0.52, 0.14, 'Energy barrier', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.60, 0.31, 'Pinning site', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.69, 0.20, 'Energy barrier', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(-0.77, 0.36, 'Pinning site', 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(0.34, 0.335, 'k_-(x+\Delta{}x)', 'HorizontalAlignment', 'Center', 'Interpreter', 'tex');
text(0.34, 0.4, 'k_+(x+\Delta{}x)', 'HorizontalAlignment', 'Center', 'Interpreter', 'tex');
text(-0.14, 0.27, 'k_-(x)', 'HorizontalAlignment', 'Center', 'Interpreter', 'tex');
text(-0.11, 0.37, 'k_+(x)', 'HorizontalAlignment', 'Center', 'Interpreter', 'tex');

arrowT = linspace(140, 240, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.05*cos(1+1.69*arrowt) + 0.05*arrowt + 0.12; 
aE = arrowE(1); 
aT = (arrowT(1)-140)/800+0.1; 

plot((arrowT-140)/800+0.1, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT+0.04 aT], [aE aE+0.02 aE+0.03], 'k'); 


arrowT = linspace(140, 240, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.05*cos(1+1.69*arrowt) + 0.05*arrowt + 0.17; 
aE = arrowE(end); 
aT = (arrowT(end)-140)/800+0.1; 

plot((arrowT-140)/800+0.1, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT-0.025 aT-0.05], [aE aE+0.025 aE-0.0], 'k'); 


arrowT = linspace(140, 240, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.05*cos(1+1.69*arrowt) + 0.05*arrowt + 0.09; 
aE = arrowE(1); 
aT = (arrowT(1)-140)/800-0.08; 

plot((arrowT-140)/800-0.08, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT+0.04 aT], [aE aE+0.02 aE+0.03], 'k'); 


arrowT = linspace(140, 240, 1000); 
arrowt = arrowT / 180 * pi; 
arrowE = 0.05*cos(1+1.69*arrowt) + 0.05*arrowt + 0.15; 
aE = arrowE(end); 
aT = (arrowT(end)-140)/800-0.08; 

plot((arrowT-140)/800-0.08, arrowE, 'k-', 'LineWidth', 2);
patch([aT aT-0.025 aT-0.05], [aE aE+0.025 aE-0.0], 'k'); 


ax = gca;
ax.XTick = -1:0.5:1;
ax.YTick = [];


title('c) Proposed new MD theory'); 



print(gcf, '-dpng', '../output/png/Schema.png'); 
print(gcf, '-dpdf', '../output/pdf/Schema.pdf', '-bestfit'); 



