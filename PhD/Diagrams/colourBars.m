%% Colour Bars

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

levels = 32;

colormap(viridis(levels));
% colormap(flipud(viridis(levels)));
% colormap(cool2warm(levels));

% Figure Setup
fig = fig + 1;
figName = 'Colour_Bar_QS_Base_Cp';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Figure Formatting
axis off;
box off;
caxis([-0.245; -0.09]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 5):max(caxis);
cB = colorbar('ticks', tickData(2:end-1), 'location', 'southOutside', 'axisLocation', 'out');
% cB.Ruler.TickLabelFormat = '%.0f';
% cB.Ruler.TickLabelFormat = '%.1f';
% cB.Ruler.TickLabelFormat = '%.2f';
% cB.Ruler.TickLabelFormat = '%.3f';
% cB.Ruler.TickLabelFormat = '%.4f';
% cB.Ruler.TickLabelFormat = '%+.0f';
% cB.Ruler.TickLabelFormat = '%+.1f';
% cB.Ruler.TickLabelFormat = '%+.2f';
cB.Ruler.TickLabelFormat = '%+.3f';
% cB.Ruler.TickLabelFormat = '%+.4f';
cB.Ruler.Exponent = 0;
cB.Label.Interpreter = 'latex';
cB.Label.String = '{$C_{p}$}';

pause(0.5);
hold off;

% Save Figure
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');