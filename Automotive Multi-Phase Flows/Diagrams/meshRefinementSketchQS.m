%% Quarter-Scale Mesh Refinement Sketch v1.0
% ----
% Schematic Diagram of Quarter-Scale Mesh Refinement Regions


%% Preamble

run preamble;

disp('=========================================');
disp('Quarter-Scale Mesh Refinement Sketch v1.0');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Geometries

% Select Subject Geometry
[geometry, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

% Select CAD-Based Refinement Regions
[CADrefinement, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

% Define Cube-Based Reginement Regions
cubeRefinement.underbodyInnerA = [
                                  -1.92575, -0.2345, 0
                                  -1.79175 -0.2345, 0
                                  -1.79175, 0.2345, 0
                                  -1.92575, 0.2345, 0
                                  -1.92575, -0.2345, 0.125
                                  -1.79175 -0.2345, 0.125
                                  -1.79175, 0.2345, 0.125
                                  -1.92575, 0.2345, 0.125
                                  ];

cubeRefinement.underbodyInnerB = [
                                  -1.79175, -0.2345, 0
                                  -0.80175 -0.2345, 0
                                  -0.80175, 0.2345, 0
                                  -1.79175, 0.2345, 0
                                  -1.79175, -0.2345, 0.18
                                  -0.80175 -0.2345, 0.18
                                  -0.80175, 0.2345, 0.18
                                  -1.79175, 0.2345, 0.18
                                 ];

cubeRefinement.underbodyOuterA = [
                                  -1.96575, -0.2745, 0
                                  -1.79175 -0.2745, 0
                                  -1.79175, 0.2745, 0
                                  -1.96575, 0.2745, 0
                                  -1.96575, -0.2745, 0.125
                                  -1.79175 -0.2745, 0.125
                                  -1.79175, 0.2745, 0.125
                                  -1.96575, 0.2745, 0.125
                                 ];

cubeRefinement.underbodyOuterB = [
                                  -1.79175, -0.2745, 0
                                  -0.864 -0.2745, 0
                                  -0.864, 0.2745, 0
                                  -1.79175, 0.2745, 0
                                  -1.79175, -0.2745, 0.18
                                  -0.864 -0.2745, 0.18
                                  -0.864, 0.2745, 0.18
                                  -1.79175, 0.2745, 0.18
                                 ];

cubeRefinement.wheelsFL = [
                           -1.7, -0.2065, 0
                           -1.5875 -0.2065, 0
                           -1.5875, -0.1275, 0
                           -1.7, -0.1275, 0
                           -1.7, -0.2065, 0.016
                           -1.5875 -0.2065, 0.016
                           -1.5875, -0.1275, 0.016
                           -1.7, -0.1275, 0.016
                         ];

cubeRefinement.wheelsFR = [
                           -1.7, 0.1275, 0
                           -1.5875 0.1275, 0
                           -1.5875, 0.2065, 0
                           -1.7, 0.2065, 0
                           -1.7, 0.1275, 0.016
                           -1.5875 0.1275, 0.016
                           -1.5875, 0.2065, 0.016
                           -1.7, 0.2065, 0.016
                         ];

cubeRefinement.wheelsRL = [
                           -1.0625, -0.2065, 0
                           -0.95 -0.2065, 0
                           -0.95, -0.1275, 0
                           -1.0625, -0.1275, 0
                           -1.0625, -0.2065, 0.016
                           -0.95 -0.2065, 0.016
                           -0.95, -0.1275, 0.016
                           -1.0625, -0.1275, 0.016
                         ];

cubeRefinement.wheelsRR = [
                           -1.0625, 0.1275, 0
                           -0.95 0.1275, 0
                           -0.95, 0.2065, 0
                           -1.0625, 0.2065, 0
                           -1.0625, 0.1275, 0.016
                           -0.95 0.1275, 0.016
                           -0.95, 0.2065, 0.016
                           -1.0625, 0.2065, 0.016
                         ];

cubeConnect = [
               4, 8, 5, 1, 4
               1, 5, 6, 2, 1
               2, 6, 7, 3, 2
               3, 7, 8, 4, 3
               5, 8, 7, 6, 5
               1, 4, 3, 2, 1
              ]';

% Define Domain Vertices
domainVerts = [
                -9.44, -0.9519083, 0 % 0
                -4.736, -0.9519083, 0 % 1
                -4.608, -0.9519083, 0 % 2
                4.832, -0.9781306, 0 % 3
                4.960, -0.9781306, 0 % 4
                9.664, -0.9781306, 0 % 5
                -9.44, -0.9519083, 1.32 % 6
                -4.736, -0.9519083, 1.32 % 7
                -4.608, -0.9519083, 1.32 % 8
                4.832, -0.9781306, 1.32 % 9
                4.960, -0.9781306, 1.32 % 10
                9.664, -0.9781306, 1.32 % 11
                -9.44, 0.9519083, 0 % 12
                -4.736, 0.9519083, 0 % 13
                -4.608, 0.9519083, 0 % 14
                4.832, 0.9781306, 0 % 15
                4.960, 0.9781306, 0 % 6
                9.664, 0.9781306, 0 % 17
                -9.44, 0.9519083, 1.32 % 18
                -4.736, 0.9519083, 1.32 % 19
                -4.608, 0.9519083, 1.32 % 20
                4.832, 0.9781306, 1.32 % 21
                4.960, 0.9781306, 1.32 % 22
                9.664, 0.9781306, 1.32 % 23
               ];

% Define Vertex Connectivity
wsConnect = [
             14, 15, 21, 20, 14 % Back
             2, 3, 15, 14, 2 % Floor
             2, 3, 9, 8, 2 % Front
             8, 9, 21, 20, 8 % Roof
            ] + 1;


%% Present Refinement Regions

% Initialise Figure
fig = fig + 1;
figName = 'Quarter_Scale_Mesh_Refinement_Regions';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
lighting gouraud;

% Plot Geometry
if ~isempty(geometry)
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        patch('faces', geometry.(parts{i}).faces, ...
              'vertices', geometry.(parts{i}).vertices, ...
              'faceColor', [0.5, 0.5, 0.5], ...
              'edgeColor', [0.5, 0.5, 0.5], ...
              'lineStyle', 'none');
    end
    clear i parts;

end

% Plot CAD-Based Refinement Regions
if ~isempty(CADrefinement)
    
    parts = fieldnames(CADrefinement);
    for i = 1:height(parts)
        patch('faces', CADrefinement.(parts{i}).faces, ...
              'vertices', CADrefinement.(parts{i}).vertices, ...
              'faceColor', graphColours(1), ...
              'edgeColor', graphColours(1), ...
              'faceAlpha', 0.25, ...
              'lineStyle', 'none');
    end
    clear i parts opacity;

end

% Plot Cube-Based Reginement Regions
regions = fieldnames(cubeRefinement);
for i = 1:height(regions)
    x = cubeRefinement.(regions{i})(:,1);
    y = cubeRefinement.(regions{i})(:,2);
    z = cubeRefinement.(regions{i})(:,3);
    
    patch('xData', x(cubeConnect), 'yData', y(cubeConnect), 'zData', z(cubeConnect), ...
                                                            'faceColor', graphColours(1), ...
                                                            'edgeColor', graphColours(1), ...
                                                            'faceAlpha', 0.25, ...
                                                            'lineStyle', 'none');
end
clear i x y z regions;

% Draw Domain
for i = 1:height(wsConnect)
    plot3(domainVerts(wsConnect(i,:),1), ...
          domainVerts(wsConnect(i,:),2), ...
          domainVerts(wsConnect(i,:),3), ...
          'color', graphColours(2), 'lineWidth', 2);
end
clear i;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }', 'interpreter', 'latex');
axis off;
view([-60, 15]);
lightangle(-60,15)
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');