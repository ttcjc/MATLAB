%% Full-Scale Computational Domain Sketch v1.0
% ----
% Schematic Diagram of Computational Domain Representative of On-Road Conditions


%% Preamble

run preamble;

disp('===========================================');
disp('Full-Scale Computational Domain Sketch v1.0');
disp('===========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Geometries

% Select Subject Geometry
[geometry, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

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


%% Present Domain Sketch

% Initialise Figure
fig = fig + 1;
figName = 'Full_Scale_Computational_Domain';
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
    clear i;

end

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