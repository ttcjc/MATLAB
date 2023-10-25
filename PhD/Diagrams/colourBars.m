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
% colormap(jet(levels));

% Figure Setup
fig = fig + 1;
figName = 'Colour_Bar_Isolated_Wheel_U';
set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Figure Formatting
axis off;
box off;
caxis([0; 25]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 5):max(caxis);
cB = colorbar('ticks', tickData(2:end-1), 'location', 'south', 'axisLocation', 'out');
% cB.Ruler.TickLabelFormat = '%.0f';
% cB.Ruler.TickLabelFormat = '%.1f';
% cB.Ruler.TickLabelFormat = '%.2f';
% cB.Ruler.TickLabelFormat = '%.3f';
% cB.Ruler.TickLabelFormat = '%.4f';
% cB.Ruler.TickLabelFormat = '%+.0f';
% cB.Ruler.TickLabelFormat = '%+.1f';
% cB.Ruler.TickLabelFormat = '%+.2f';
% cB.Ruler.TickLabelFormat = '%+.3f';
cB.Ruler.Exponent = 0;
cB.Label.String = {'Velocity Magnitude'; '[\it{m\cdots^{-1}}]'};
% cB.Label.String = {'Velocity Magnitude'; '[{\it{m\cdots^{-1}}}]'};
hold off;

% Save Plot
pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);