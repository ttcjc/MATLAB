close all;
clc;


%%

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;


%%

xLimsPlot = [0.3; 4.6257662];
yLimsPlot = [-0.5; 0.5];
zLimsPlot = [0; 0.5];

positionData = LagData.positionCartesian{100} / normLength;

figTime = num2str(LagData.time(100), ['%.', num2str(timePrecision), 'f']);
figTitle = ['{', figTime, ' \it{s}}'];


%%

fig = fig + 1;
figName = 'Discrete_Parcel_Locations';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
lighting gouraud;

% Plot Geometry
parts = fieldnames(geometry);
for i = 1:height(parts)
    patch('faces', geometry.(parts{i}).faces, ...
          'vertices', geometry.(parts{i}).vertices, ...
          'faceColor', [0.5, 0.5, 0.5], ...
          'edgeColor', [0.5, 0.5, 0.5], ...
          'lineStyle', 'none');
end
clear i parts;

% Plot Parcels
plot3(positionData(:,1), positionData(:,2), positionData(:,3), ...
      'lineStyle', 'none', 'marker', 'o', 'markerSize', 1, 'markerFaceColor', graphColours(1), ...
      'markerEdgeColor', graphColours(1));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle(figTitle);
lightangle(90, 45);
axis on;
box on;
grid off;
view([90, 0]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = yLimsPlot(1):(diff(yLimsPlot) / 5):yLimsPlot(2);
yticks(tickData(2:5));
tickData = zLimsPlot(1):(diff(zLimsPlot) / 5):zLimsPlot(2);
zticks(tickData(2:5));
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
ylabel({'{$y_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
zlabel({'{-----}'; '{$z_{_{\ell}}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
