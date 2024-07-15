%% Full-Scale Mesh Refinement Sketch v1.0
% ----
% Schematic Diagram of Full-Scale Mesh Refinement Regions


%% Preamble

run preamble;

disp('======================================');
disp('Full-Scale Mesh Refinement Sketch v1.0');
disp('======================================');

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
                                  -2.403, -0.938, 0
                                  -1.867, -0.938, 0
                                  -1.867, 0.938, 0
                                  -2.403, 0.938, 0
                                  -2.403, -0.938, 0.482
                                  -1.867, -0.938, 0.482
                                  -1.867, 0.938, 0.482
                                  -2.403, 0.938, 0.482
                                  ];

cubeRefinement.underbodyInnerB = [
                                  -1.867, -0.938, 0
                                  2.093, -0.938, 0
                                  2.093, 0.938, 0
                                  -1.867, 0.938, 0
                                  -1.867, -0.938, 0.702
                                  2.093, -0.938, 0.702
                                  2.093, 0.938, 0.702
                                  -1.867, 0.938, 0.702
                                 ];

cubeRefinement.underbodyOuterA = [
                                  -1.96575, -1.098, 0
                                  -1.867, -1.098, 0
                                  -1.867, 1.098, 0
                                  -1.96575, 1.098, 0
                                  -1.96575, -1.098, 0.482
                                  -1.867, -1.098, 0.482
                                  -1.867, 1.098, 0.482
                                  -1.96575, 1.098, 0.482
                                 ];

cubeRefinement.underbodyOuterB = [
                                  -1.867, -1.098, 0
                                  1.774, -1.098, 0
                                  1.774, 1.098, 0
                                  -1.867, 1.098, 0
                                  -1.867, -1.098, 0.702
                                  1.774, -1.098, 0.702
                                  1.774, 1.098, 0.702
                                  -1.867, 1.098, 0.702
                                 ];

cubeRefinement.wheelsFL = [
                           -1.5, -0.826, 0
                           -1.05, -0.826, 0
                           -1.05, -0.51, 0
                           -1.5, -0.51, 0
                           -1.5, -0.826, 0.046
                           -1.05, -0.826, 0.046
                           -1.05, -0.51, 0.046
                           -1.5, -0.51, 0.046
                         ];

cubeRefinement.wheelsFR = [
                           -1.5, 0.51, 0
                           -1.05, 0.51, 0
                           -1.05, 0.826, 0
                           -1.5, 0.826, 0
                           -1.5, 0.51, 0.046
                           -1.05, 0.51, 0.046
                           -1.05, 0.826, 0.046
                           -1.5, 0.826, 0.046
                         ];

cubeRefinement.wheelsRL = [
                           1.05, -0.826, 0
                           1.5, -0.826, 0
                           1.5, -0.51, 0
                           1.05, -0.51, 0
                           1.05, -0.826, 0.046
                           1.5, -0.826, 0.046
                           1.5, -0.51, 0.046
                           1.05, -0.51, 0.046
                         ];

cubeRefinement.wheelsRR = [
                           1.05, 0.51, 0
                           1.5, 0.51, 0
                           1.5, 0.826, 0
                           1.05, 0.826, 0
                           1.05, 0.51, 0.046
                           1.5, 0.51, 0.046
                           1.5, 0.826, 0.046
                           1.05, 0.826, 0.046
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
               -27.648, -13.312, 0 % 0
               60.416, -13.312, 0 % 1
               -27.648, -13.312, 14.336 % 2
               60.416, -13.312, 14.336 % 3
               -27.648, 13.312, 0 % 4
               60.416, 13.312, 0 % 5
               -27.648, 13.312, 14.336 % 6
               60.416, 13.312, 14.336 % 7
              ];

% Define Vertex Connectivity
domainConnect = [
                 4, 5, 7, 6, 4 % Back
                 0, 1, 5, 4, 0 % Floor
                 0, 1, 3, 2, 0 % Front
                 2, 3, 7, 6, 2 % Roof
                ] + 1;


%% Present Refinement Regions

% Initialise Figure
fig = fig + 1;
figName = 'Full_Scale_Mesh_Refinement_Regions';
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
for i = 1:height(domainConnect)
    plot3(domainVerts(domainConnect(i,:),1), ...
          domainVerts(domainConnect(i,:),2), ...
          domainVerts(domainConnect(i,:),3), ...
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