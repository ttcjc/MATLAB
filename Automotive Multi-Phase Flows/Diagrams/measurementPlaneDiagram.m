clear variables;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(false);

close all
clear measurementPlanes

% Define Wind Tunnel Bounding Box
LWT.vertices = [
                -3, -0.81, 0;
                -3, 0.81, 0;
                -3, 0.81, 1.32;
                -3, -0.81, 1.32;
                3, -0.81, 0;
                3, 0.81, 0;
                3, 0.81, 1.32;
                3, -0.81, 1.32;
           ];

LWT.connectivity = [
                    1, 5, 8, 4, 1; % Front
                    2, 6, 7, 3, 2; % Back
                    4, 8, 7, 3, 4; % Roof
                    1, 5, 6, 2, 1; % Floor
                   ];

% Define Measurement Plane(s)
measurementPlanes.A.position = [
                                0.20225, -0.43, 0.01;
                                0.20225, -0.43, 0.4995;
                                0.20225, 0.25, 0.4995;
                                0.20225, 0.25, 0.01;
                               ];
measurementPlanes.A.name = '1.0 $\ell$';

measurementPlanes.B.position = [
                                0.72425, -0.43, 0.01;
                                0.72425, -0.43, 0.4995;
                                0.72425, 0.25, 0.4995;
                                0.72425, 0.25, 0.01;
                               ];
measurementPlanes.B.name = '1.5 $\ell$';

measurementPlanes.C.position = [
                                1.24625, -0.43, 0.01;
                                1.24625, -0.43, 0.4995;
                                1.24625, 0.25, 0.4995;
                                1.24625, 0.25, 0.01;
                               ];
measurementPlanes.C.name = '2.0 $\ell$';

% % Define Measurement Plane(s)
% measurementPlanes.A.position = [
%                                 0.20225, -0.522, 0;
%                                 0.20225, -0.522, 0.6264;
%                                 0.20225, 0.522, 0.6264;
%                                 0.20225, 0.522, 0;
%                                ];
% measurementPlanes.A.name = '1.0 $\ell$';
% 
% measurementPlanes.B.position = [
%                                 0.72425, -0.522, 0;
%                                 0.72425, -0.522, 0.6264;
%                                 0.72425, 0.522, 0.6264;
%                                 0.72425, 0.522, 0;
%                                ];
% measurementPlanes.B.name = '1.5 $\ell$';
% 
% measurementPlanes.C.position = [
%                                 1.24625, -0.522, 0;
%                                 1.24625, -0.522, 0.6264;
%                                 1.24625, 0.522, 0.6264;
%                                 1.24625, 0.522, 0;
%                                ];
% measurementPlanes.C.name = '2.0 $\ell$';

% Figure Settings
figName = 'Methodology_Measurement_Planes';
cMap = viridis(height(fieldnames(measurementPlanes)));
% cMap = viridis(3);
% cMap = cMap(1,:);
xLimsPlot = [-1.00625, 1.41075];
yLimsPlot = [-0.81, 0.81];
zLimsPlot = [0, 0.81];

% Figure Setup
fig = fig + 1;
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
lighting gouraud;
hold on;

% Figure Plotting
parts = fieldnames(geometry);
planes = fieldnames(measurementPlanes);

for i = 1:height(parts)
    patch('faces', geometry.(parts{i}).faces, ...
          'vertices', geometry.(parts{i}).vertices, ...
          'faceColor', ([128, 128, 128] / 255), ...
          'edgeColor', ([128, 128, 128] / 255), ...
          'lineStyle', 'none');
end

for i = 1:height(planes)
    patch(measurementPlanes.(planes{i}).position(:,1), ...
          measurementPlanes.(planes{i}).position(:,2), ...
          measurementPlanes.(planes{i}).position(:,3), ...
          cMap(i,:), ...
          'lineWidth', 2, ...
          'faceLighting', 'none');
    
    text(measurementPlanes.(planes{i}).position(3,1), ...
         measurementPlanes.(planes{i}).position(3,2), ...
         measurementPlanes.(planes{i}).position(3,3), ...
         measurementPlanes.(planes{i}).name, ...
         'interpreter', 'latex', ...
         'fontName', 'LM Mono 12', ...
         'fontSize', 20, ...
         'fontWeight', 'bold', ...
         'rotation', -90, ...
         'horizontalAlignment', 'left', ...
         'verticalAlignment', 'bottom')
end

for i = 1:height(LWT.connectivity)
    plot3(LWT.vertices(LWT.connectivity(i,:),1), ...
          LWT.vertices(LWT.connectivity(i,:),2), ...
          LWT.vertices(LWT.connectivity(i,:),3), ...
          'color', 'k', ...
          'lineWidth', 4)
end

% Figure Formatting
lightangle(0, 45);
axis off;
box off;
view([25, 15])
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = [];
yticks(tickData);
tickData = [];
zticks(tickData);
hold off;

pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);