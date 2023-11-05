clc;
close all;

% Figure Setup
fig = fig + 1;
set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', 'Rays');
set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
lighting gouraud;
hold on;

% Geometry Plotting
if ~isempty(geometry)
    parts = fieldnames(geometry);

    for i = 1:height(parts)
        patch('faces', geometry.(parts{i}).faces, ...
              'vertices', geometry.(parts{i}).vertices, ...
              'faceColor', [0.5, 0.5, 0.5], ...
              'edgeColor', [0.5, 0.5, 0.5], ...
              'lineStyle', 'none');
    end

end

y = unique(LOSdata.positionGrid(:,2));
y = y(20:20:140);
z = unique(LOSdata.positionGrid(:,3));
z = z(20:20:80);
x = ones(size(z)) * LOSdata.positionGrid(1,1);

[x,y,z] = ndgrid(x,y,z);

positions = [x(:), y(:), z(:)];

plane = [
         LOSdata.positionGrid(1,1), min(LOSdata.positionGrid(:,2)), min(LOSdata.positionGrid(:,3));
         LOSdata.positionGrid(1,1), min(LOSdata.positionGrid(:,2)), max(LOSdata.positionGrid(:,3));
         LOSdata.positionGrid(1,1), max(LOSdata.positionGrid(:,2)), max(LOSdata.positionGrid(:,3));
         LOSdata.positionGrid(1,1), max(LOSdata.positionGrid(:,2)), min(LOSdata.positionGrid(:,3));
         LOSdata.positionGrid(1,1), min(LOSdata.positionGrid(:,2)), min(LOSdata.positionGrid(:,3));
        ];

LOSdata.originPoint(2) = 0;

% Data Plotting
scatter3(LOSdata.originPoint(1), LOSdata.originPoint(2), LOSdata.originPoint(3), 100, ([230, 0, 126] / 255), 'filled');
plot3(plane(:,1), plane(:,2), plane(:,3), 'color', ([230, 0, 126] / 255), 'lineWidth', 2);

for i = 1:height(positions)
    plot3([LOSdata.originPoint(1); positions(i,1)], ...
          [LOSdata.originPoint(2); positions(i,2)], ...
          [LOSdata.originPoint(3); positions(i,3)], ...
          'color', ([74, 24, 99, 127.5] / 255), 'lineWidth', 2);
end

% Figure Formatting
title(figTitle, 'color', ([254, 254, 254] / 255));
subtitle(figSubtitle);
lightangle(0, 45);
axis on;
box on;
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
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);
exportgraphics(gca, '~/MATLAB/Output/Figures/Rays.png', 'resolution', 300);
savefig(gcf, '~/MATLAB/Output/Figures/Rays.fig')