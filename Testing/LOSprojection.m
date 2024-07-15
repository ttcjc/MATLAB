run preamble;

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%%

xLimsData = [-0.537116858237548; 4.462883141762452];
yLimsData = [-0.6; 0.6];
zLimsData = [0; 0.6];

xLimsPlot = [-0.637116858237548; 4.562883141762452];
yLimsPlot = [-0.7; 0.7];
zLimsPlot = [0; 0.7];

planePosition = 1.949 / normLength;

positionsPlane = [
                  planePosition, yLimsData(1), zLimsData(1);
                  planePosition, yLimsData(1), zLimsData(2);
                  planePosition, yLimsData(2), zLimsData(2);
                  planePosition, yLimsData(2), zLimsData(1);
                  planePosition, yLimsData(1), zLimsData(1);
                 ];

positionsSensor = [18.637, 0, 0.76] / normLength;

positionsRay = dsearchn((volumeData.positionGrid(:,[2,3]) / normLength), ([0, 0.76] / normLength));
positionsRay = [planePosition, (volumeData.positionGrid(positionsRay, [2,3]) / normLength)];


%%

% Initialise Figure
fig = fig + 1;
figName = 'LoS_Projection';
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
          'vertices', (geometry.(parts{i}).vertices / normLength), ...
          'faceColor', [0.5, 0.5, 0.5], ...
          'edgeColor', [0.5, 0.5, 0.5], ...
          'lineStyle', 'none');
end
clear i parts;

% Plot Ray Projection
plot3(positionsSensor(:,1), positionsSensor(:,2), positionsSensor(:,3), 'lineStyle', 'none', ...
                                                                        'marker', 'o', ...
                                                                        'markerSize', 10, ...
                                                                        'markerFaceColor', graphColours(2), ...
                                                                        'markerEdgeColor', graphColours(2));

plot3(positionsPlane(:,1), positionsPlane(:,2), positionsPlane(:,3), 'lineStyle', '-', ...
                                                                     'lineWidth', 2, ...
                                                                     'color', graphColours(2));

for i = 1:height(positionsRay)
    plot3([positionsSensor(:,1); positionsRay(i,1)], ...
          [positionsSensor(:,2); positionsRay(i,2)], ... 
          [positionsSensor(:,3); positionsRay(i,3)], 'lineStyle', '-', ...
                                                     'lineWidth', 4, ...
                                                     'color', ([74, 24, 99] / 255));
end
clear i;

% for i = 1:4
%     plot3([positionsSensor(:,1); positionsPlane(i,1)], ...
%           [positionsSensor(:,2); positionsPlane(i,2)], ... 
%           [positionsSensor(:,3); positionsPlane(i,3)], 'lineStyle', '-', ...
%                                                        'lineWidth', 2, ...
%                                                        'color', ([74, 24, 99, 127.5] / 255));
% end
% clear i;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
lightangle(0, 45);
axis on;
box on;
grid off;
view([30, 30]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = [];
yticks(tickData);
tickData = [];
zticks(tickData);
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');