run preamble;

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

% Define Wind Tunnel Bounding Box
LWT.vertices = [
                -3, -0.81, 0;
                -3, 0.81, 0;
                -3, 0.81, 1.32;
                -3, -0.81, 1.32;
                4, -0.81, 0;
                4, 0.81, 0;
                4, 0.81, 1.32;
                4, -0.81, 1.32;
           ];

LWT.connectivity = [
                    1, 5, 8, 4, 1; % Front
                    2, 6, 7, 3, 2; % Back
                    4, 8, 7, 3, 4; % Roof
                    1, 5, 6, 2, 1; % Floor
                   ];

% % Define Experimental Measurement Plane(s)
measurementPlanes.A.position = [
                                0.20225, -0.43, 0.01;
                                0.20225, -0.43, 0.4995;
                                0.20225, 0.25, 0.4995;
                                0.20225, 0.25, 0.01;
                               ];
measurementPlanes.A.name = '1.0 $\ell$';
measurementPlanes.A.colour = graphColours(1);

measurementPlanes.B.position = [
                                0.72425, -0.43, 0.01;
                                0.72425, -0.43, 0.4995;
                                0.72425, 0.25, 0.4995;
                                0.72425, 0.25, 0.01;
                               ];
measurementPlanes.B.name = '1.5 $\ell$';
measurementPlanes.B.colour = graphColours(2);

measurementPlanes.C.position = [
                                1.24625, -0.43, 0.01;
                                1.24625, -0.43, 0.4995;
                                1.24625, 0.25, 0.4995;
                                1.24625, 0.25, 0.01;
                               ];
measurementPlanes.C.name = '2.0 $\ell$';
measurementPlanes.C.colour = graphColours(3);


% Define Numerical Measurement Plane(s)
% measurementPlanes.A.position = [
%                                 0.20225, -0.522, 0;
%                                 0.20225, -0.522, 0.522;
%                                 0.20225, 0.522, 0.522;
%                                 0.20225, 0.522, 0;
%                                ];
% measurementPlanes.A.name = '1.0 $\ell$';
% measurementPlanes.A.colour = graphColours(1);

% measurementPlanes.B.position = [
%                                 0.72425, -0.522, 0;
%                                 0.72425, -0.522, 0.522;
%                                 0.72425, 0.522, 0.522;
%                                 0.72425, 0.522, 0;
%                                ];
% measurementPlanes.B.name = '1.5 $\ell$';
% measurementPlanes.B.colour = graphColours(2);

% measurementPlanes.C.position = [
%                                 1.24625, -0.522, 0;
%                                 1.24625, -0.522, 0.522;
%                                 1.24625, 0.522, 0.522;
%                                 1.24625, 0.522, 0;
%                                ];
% measurementPlanes.C.name = '2.0 $\ell$';
% measurementPlanes.C.colour = graphColours(3);
% 
% measurementPlanes.D.position = [
%                                 2.29025, -0.522, 0;
%                                 2.29025, -0.522, 0.522;
%                                 2.29025, 0.522, 0.522;
%                                 2.29025, 0.522, 0;
%                                ];
% measurementPlanes.D.name = '3.0 $\ell$';
% measurementPlanes.D.colour = graphColours(4);

% measurementPlanes.E.position = [
%                                 3.33425, -0.522, 0;
%                                 3.33425, -0.522, 0.522;
%                                 3.33425, 0.522, 0.522;
%                                 3.33425, 0.522, 0;
%                                ];
% measurementPlanes.E.name = '4.0 $\ell$';
% measurementPlanes.E.colour = graphColours(5);

% Figure Settings
figName = 'Methodology_Measurement_Planes';
% cMap = viridis(height(fieldnames(measurementPlanes)));
% cMap = viridis(3);
% cMap = cMap(1,:);
xLimsPlot = [-1.00625, 1.41075]; % xLimsPlot = [-1.00625, 3.49875];
yLimsPlot = [-0.81, 0.81];
zLimsPlot = [0, 0.81];

% Figure Setup
fig = fig + 1;
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880])
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
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
    patch('xData', measurementPlanes.(planes{i}).position(:,1), ...
          'yData', measurementPlanes.(planes{i}).position(:,2), ...
          'zData', measurementPlanes.(planes{i}).position(:,3), ...
          'faceColor', measurementPlanes.(planes{i}).colour, ...
          'lineWidth', 2, ...
          'faceLighting', 'none');
    
    text(measurementPlanes.(planes{i}).position(3,1), ...
         measurementPlanes.(planes{i}).position(3,2), ...
         measurementPlanes.(planes{i}).position(3,3), ...
         measurementPlanes.(planes{i}).name, ...
         'interpreter', 'latex', ...
         'fontName', 'LM Mono 12', ...
         'fontSize', 22, ...
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
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');