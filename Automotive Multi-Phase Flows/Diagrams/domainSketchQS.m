%% Quarter-Scale Computational Domain Sketch v1.0
% ----
% Schematic Diagram of Computational Domain Representative of the Loughborough University


%% Preamble

run preamble;

disp('==============================================');
disp('Quarter-Scale Computational Domain Sketch v1.0');
disp('==============================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Geometries

% Select Subject Geometry
[geometry, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

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
inletConnect = [
                12, 13, 19, 18, 12 % Back
                0, 1, 13, 12, 0 % Floor
                0, 1, 7, 6, 0 % Front
                6, 7, 19, 18, 6 % Roof
               ] + 1;

inletBufferConnect = [
                      13, 14, 20, 19, 13 % Back
                      1, 2, 14, 13, 1 % Floor
                      1, 2, 8, 7, 1 % Front
                      7, 8, 20, 19, 7 % Roof
                     ] + 1;

workingSectionConnect = [
                         14, 15, 21, 20, 14 % Back
                         2, 3, 15, 14, 2 % Floor
                         2, 3, 9, 8, 2 % Front
                         8, 9, 21, 20, 8 % Roof
                        ] + 1;

outletBufferConnect = [
                       15, 16, 22, 21, 15 % Back
                       3, 4, 16, 15, 3 % Floor
                       3, 4, 10, 9, 3 % Front
                       9, 10, 22, 21, 9 % Roof
                      ] + 1;

outletConnect = [
                 16, 17, 23, 22, 16 % Back
                 4, 5, 17, 16, 4 % Floor
                 4, 5, 11, 10, 4 % Front
                 10, 11, 23, 22, 10 % Roof
                ] + 1;


%% Present Domain Sketch

% Initialise Figure
fig = fig + 1;
figName = 'Quarter_Scale_Computational_Domain';
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
regions = ["inletConnect"; "workingSectionConnect"; "outletConnect"; "inletBufferConnect"; "outletBufferConnect"];

legEntries = zeros([height(regions), 1]);

colours = [graphColours(1); graphColours(2); graphColours(3); graphColours(7); graphColours(7)];
for i = 1:height(regions)
    region = eval(regions(i));
    
    for j = 1:height(region)
        
        if j == 1
            legEntries(i) = plot3(domainVerts(region(j,:),1), ...
                                  domainVerts(region(j,:),2), ...
                                  domainVerts(region(j,:),3), ...
                                  'color', colours(i,:), 'lineWidth', 2);
        else
            plot3(domainVerts(region(j,:),1), ...
                  domainVerts(region(j,:),2), ...
                  domainVerts(region(j,:),3), ...
                  'color', colours(i,:), 'lineWidth', 2);
        end
        
    end
    
end
clear colours;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }', 'interpreter', 'latex');
axis off;
view([-60, 15]);
lightangle(-60,15)
legend(legEntries([1, 2, 3, 4]), 'Inlet', 'Working Section', 'Outlet', 'Buffer Region', ...
       'position', [-0.2415, 0.137, 1, 1], 'interpreter', 'latex');
legend boxoff;
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');