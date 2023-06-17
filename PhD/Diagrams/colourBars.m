%% Colour Bars

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

% colormap(viridis(32));
% colormap(flipud(viridis(32)));
colormap(cool2warm(32));

% Figure Setup
fig = fig + 1;
figName = 'Colour_Bar_POD';
set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Figure Formatting
axis off;
box off;
caxis([-1; 1]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 5):max(caxis);
cB = colorbar('ticks', tickData(2:end-1), 'location', 'east', 'axisLocation', 'out');
% cB.Ruler.TickLabelFormat = '%.0f';
% cB.Ruler.TickLabelFormat = '%.1f';
% cB.Ruler.TickLabelFormat = '%.2f';
% cB.Ruler.TickLabelFormat = '%.3f';
% cB.Ruler.TickLabelFormat = '%.4f';
% cB.Ruler.TickLabelFormat = '%+.0f';
cB.Ruler.TickLabelFormat = '%+.1f';
% cB.Ruler.TickLabelFormat = '%+.2f';
% cB.Ruler.TickLabelFormat = '%+.3f';
cB.Ruler.Exponent = 0;
cB.Label.String = '\Phi_n';
% cB.Label.String = {'Contaminant Mass Area Density'; '[{\it{g\cdotm^{-2}}}]'};
hold off;

% Save Plot
pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);