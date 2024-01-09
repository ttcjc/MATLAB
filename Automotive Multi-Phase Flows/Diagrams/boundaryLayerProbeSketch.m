%% Boundary Layer Probe Sketch v1.0
% ----
% Lorem Ipsum


%% Preamble

run preamble;

disp('================================');
disp('Boundary Layer Probe Sketch v1.0');
disp('================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Geometries

% Select Subject Geometry
[geometry, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

% Define Domain Vertices
vertices = [
            -9.44, -0.9519083, 0; % 0
            -4.736, -0.9519083, 0; % 1
            -4.608, -0.9519083, 0; % 2
            4.832, -0.9781306, 0; % 3
            4.960, -0.9781306, 0; % 4
            9.664, -0.9781306, 0; % 5
            -9.44, -0.9519083, 1.32; % 6
            -4.736, -0.9519083, 1.32; % 7
            -4.608, -0.9519083, 1.32; % 8
            4.832, -0.9781306, 1.32; % 9
            4.960, -0.9781306, 1.32; % 10
            9.664, -0.9781306, 1.32; % 11
            -9.44, 0.9519083, 0; % 12
            -4.736, 0.9519083, 0; % 13
            -4.608, 0.9519083, 0; % 14
            4.832, 0.9781306, 0; % 15
            4.960, 0.9781306, 0; % 6
            9.664, 0.9781306, 0; % 17
            -9.44, 0.9519083, 1.32; % 18
            -4.736, 0.9519083, 1.32; % 19
            -4.608, 0.9519083, 1.32; % 20
            4.832, 0.9781306, 1.32; % 21
            4.960, 0.9781306, 1.32; % 22
            9.664, 0.9781306, 1.32 % 23
           ];

% Define Vertex Connectivity
workingSection = [
                  14, 15, 21, 20, 14; % Back
                  2, 3, 15, 14, 2; % Floor
                  2, 3, 9, 8, 2; % Front
                  8, 9, 21, 20, 8 % Roof
                 ] + 1;


%% Present Boundary Layer Probe Locations

% Initialise Figure
fig = fig + 1;
figName = 'Boundary_Layer_Probe_Locations';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
lighting gouraud;

% Plot Subject Geometry
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
for i = 1:height(workingSection)
    plot3(vertices(workingSection(i,:),1), ...
          vertices(workingSection(i,:),2), ...
          vertices(workingSection(i,:),3), ...
          'color', graphColours(2), 'lineWidth', 2);
end
clear i;

text(-4.408, 0, 0.66, 'Inlet', 'interpreter', 'latex', 'fontSize', 16,  'rotation', 90, 'horizontalAlignment', 'center');
text(4.632, 0, 0.66, 'Outlet', 'interpreter', 'latex', 'fontSize', 16,  'rotation', 270, 'horizontalAlignment', 'center');

% Draw Probe Locations
probeLocations = [-0.56075, 0, 0.66; 0, 0, 0.66];
textLocations = [-0.56075, 0, 0.66; -0.56075, -0.72, 0.66; 0, -0.72, 0.66; 0, 0, 0.66];

scatter3(probeLocations(:,1), probeLocations(:,2), probeLocations(:,3), 50, 'markerEdgeColor', graphColours(1), 'lineWidth', 2);
plot3(textLocations(1:2,1), textLocations(1:2,2), textLocations([1,2],3), 'color', graphColours(1), 'lineWidth', 2);
plot3(textLocations(3:4,1), textLocations(3:4,2), textLocations([3,4],3), 'color', graphColours(1), 'lineWidth', 2);
text(textLocations(2,1) - 0.1, textLocations(2,2) + 0, textLocations(2,3), 'Location A', 'horizontalAlignment', 'right', 'verticalAlignment', 'baseline', 'interpreter', 'latex', 'fontSize', 16);
text(textLocations(3,1) + 0.1, textLocations(3,2) + 0, textLocations(3,3), 'Location B', 'horizontalAlignment', 'left', 'verticalAlignment', 'baseline', 'interpreter', 'latex', 'fontSize', 16);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis off;
view([0, 90]);
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