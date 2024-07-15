CLA;

run preamble;

QS = load('~/MATLAB/Testing/Dispersed Phase Processing/QS_2L_nParcels.mat').totalParcels;
FS_Uncoupled = load('~/MATLAB/Testing/Dispersed Phase Processing/uncoupled_2L_nParcels.mat').totalParcels;
FS_Coupled = load('~/MATLAB/Testing/Dispersed Phase Processing/coupled_2L_nParcels.mat').totalParcels;


%%

% Initialise Figure
fig = fig + 1;
figName = 'Instantaneous_Planar_Parcel_Count';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mass Flux
plot(QS, 'color', graphColours(1), 'lineWidth', 2);
plot(FS_Uncoupled, 'color', graphColours(4), 'lineWidth', 2);
plot(FS_Coupled, 'color', graphColours(5), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 1100]);
ylim([0; 9e4]);
tickData = (220:220:880);
xticks(tickData);
tickData = (1.8e4:1.8e4:7.2e4);
yticks(tickData);
xlabel({'{Sample}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Parcel Count}'}, 'interpreter', 'latex');
legend({'\textit{QS}', ...
        '\textit{FS Uncoupled}', ...
        '\textit{FS Coupled}'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');