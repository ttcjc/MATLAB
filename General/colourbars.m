%% Colour Bars

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

% Figure Setup
fig = fig + 1;
figName = 'Colour_Bar_POD';
set(figure(fig), 'outerPosition', [25, 75, 850, 850], 'name', figName);
set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
         'fontSize', 18, 'layer', 'top');
colormap(turbo(24));
hold on;

% Figure Formatting
axis off;
box off;
caxis([-1, 1]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 5):max(caxis);
cB = colorbar('ticks', tickData(2:end-1), 'location', 'south', 'axisLocation', 'out');
cB.Ruler.TickLabelFormat = '%.1f';
cB.Ruler.Exponent = 0;
cB.Label.String = {'Normalised Mode Eigenvector', '---------'};
% cB.Label.FontName = 'LM Roman 12';
hold off;

% Save Plot
pause(2);
exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);