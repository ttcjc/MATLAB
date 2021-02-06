%% Domain Sketch

clear variables;
close all;
clc;

fig = 0;
figHold = 0; %#ok<*NASGU>

disp ('==================');
disp ('Domain Sketch v1.0');
disp ('==================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialisation

% Vertices
vertices = [-13.6, -0.96, 0;
            -4.832, -0.96, 0;
            -4.704, -0.96, 0;
            4.704, -0.98613, 0;
            4.832, -0.98613, 0;
            14.4, -0.98613, 0;
            -13.6, -0.96, 1.32;
            -4.832, -0.96, 1.32;
            -4.704, -0.96, 1.32;
            4.704, -0.98613, 1.32;
            4.832, -0.98613, 1.32;
            14.4, -0.98613, 1.32;
            -13.6, 0.96, 0;
            -4.832, 0.96, 0;
            -4.704, 0.96, 0;
            4.704, 0.98613, 0;
            4.832, 0.98613, 0;
            14.4, 0.98613, 0;
            -13.6, 0.96, 1.32;
            -4.832, 0.96, 1.32;
            -4.704, 0.96, 1.32;
            4.704, 0.98613, 1.32;
            4.832, 0.98613, 1.32;
            14.4, 0.98613, 1.32] / 1.044;

% Connectivity
inlet = [
         12, 13, 19, 18, 12; % Back
         0, 1, 13, 12, 0; % Floor
         0, 1, 7, 6, 0; % Front
         6, 7, 19, 18, 6 % Roof
        ] + 1;

workingSection = [
                  14, 15, 21, 20, 14; % Back
                  2, 3, 15, 14, 2; % Floor
                  2, 3, 9, 8, 2; % Front
                  8, 9, 21, 20, 8 % Roof
                 ] + 1;

outlet = [
          16, 17, 23, 22, 16; % Back
          4, 5, 17, 16, 4; % Floor
          4, 5, 11, 10, 4; % Front
          10, 11, 23, 22, 10 % Roof
         ] + 1;

inletBuffer = [
               13, 14, 20, 19, 13; % Back
               1, 2, 14, 13, 1; % Floor
               1, 2, 8, 7, 1; % Front
               7, 8, 20, 19, 7 % Roof
              ] + 1;

outletBuffer = [
                15, 16, 22, 21, 15; % Back
                3, 4, 16, 15, 3; % Floor
                3, 4, 10, 9, 3; % Front
                9, 10, 22, 21, 9 % Roof
               ] + 1;
           
% CAD
windsorBody = stlread('~/CAD/CFD Geometries/Windsor/Output/Presentation Quality/Windsor_Body_Square.stl');
windsorAxles = stlread('~/CAD/CFD Geometries/Windsor/Output/Presentation Quality/Windsor_Axles.stl');
windsorWheels = stlread('~/CAD/CFD Geometries/Windsor/Output/Presentation Quality/Windsor_Wheels.stl');

modelColour = [0.5, 0.5, 0.5];


%% Domain
 
% Figure Setup
fig = fig + 1;
figure('name', 'Domain');
hold on;
set(figure(fig), 'outerPosition', [1945, 25, 750, 750]);

% Draw
regions = ["inlet", "workingSection", "outlet", "inletBuffer", "outletBuffer"];

colours = [
           0.21176, 0.06667, 0.38824;
           0.71765, 0.00000, 0.38431;
           0.94902, 0.41569, 0.21961;
           0.00000, 0.00000, 0.00000;
           0.00000, 0.00000, 0.00000
          ];

for i = 1:size(regions,2)
    region = eval(regions(i));
    
    for j = 1:size(region,1)
        
        if j == 1
            legEntry(i) = plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        else
            plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        end
        
    end
    
end

% trisurf(windsorBody, 'edgeColor', modelColour, 'faceColor', modelColour);
% trisurf(windsorAxles, 'edgeColor', modelColour, 'faceColor', modelColour);
% trisurf(windsorWheels, 'edgeColor', modelColour, 'faceColor', modelColour);

% Figure Formatting
view([-75, 5]);
axis off;
legend(legEntry([1:3, 5]), 'Inlet', 'Working Section', 'Outlet', 'Buffer Region', 'position', [-0.235, 0.112, 1, 1]);
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top', ...
         'dataAspectRatio', [1, 1, 1]);
hold off;

print(fig, '~/MATLAB/Figures/Computational_Domain', '-dpng', '-r300');


%% Mesh Refinement Regions

% Figure Setup
fig = fig + 1;
figure('name', 'Refinement Regions');
hold on;
set(figure(fig), 'outerPosition', [1945, 100, 750, 750]);

% Draw
regions = "workingSection";

colours = [0.71765, 0.00000, 0.38431];

for i = 1:size(regions,2)
    region = eval(regions(i));
    
    for j = 1:size(region,1)
        
        if j == 1
            legEntry(i) = plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        else
            plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        end
        
    end
    
end

plotcube([0.126, 0.454, 0.125], [-0.593, -0.227, 0], 0.75, [0.21176, 0.06667, 0.38824]) % Underbody A
plotcube([0.982, 0.454, 0.18], [-0.467, -0.227, 0], 0.75, [0.21176, 0.06667, 0.38824]) % Underbody B
plotcube([0.878, 0.518, 0.403], [0.405, -0.259, 0], 0.5, [0.71765, 0.00000, 0.38431]) % Wake
plotcube([1.678, 0.646, 0.467], [0.405, -0.323, 0], 0.25, [0.94902, 0.41569, 0.21961]) % Local
plotcube([4.840, 0.902, 0.595], [-1.156, -0.451, 0], 0, [0 0 0]) % Far Field

trisurf(windsorBody, 'edgeColor', modelColour, 'faceColor', modelColour);
trisurf(windsorAxles, 'edgeColor', modelColour, 'faceColor', modelColour);
trisurf(windsorWheels, 'edgeColor', modelColour, 'faceColor', modelColour);

% Figure Formatting
view([-75, 5]);
axis off;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top', ...
         'dataAspectRatio', [1, 1, 1]);
hold off;

print(fig, '~/MATLAB/Figures/Refinement_Regions', '-dpng', '-r300');


%% Probe Locations
 
% Figure Setup
fig = fig + 1;
figure('name', 'Probe Locations');
hold on;
set(figure(fig), 'outerPosition', [1945, 175, 750, 750]);

% Draw
regions = "workingSection";

colours = [0.71765, 0.00000, 0.38431];

for i = 1:size(regions,2)
    region = eval(regions(i));
    
    for j = 1:size(region,1)
        
        if j == 1
            legEntry(i) = plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        else
            plot3(vertices(region(j,1:5),1), vertices(region(j,1:5),2), vertices(region(j,1:5),3), 'color', colours(i,:));
        end
        
    end
    
end

probeLocations = [-0.56075, 0, 0.66 ; 0, 0, 0.66] / 1.044;
textLocations = [
                 -0.56075, 0, 0.66 ; -0.56075, 0.72, 0.66 ;
                 0, -0.72, 0.66 ; 0, 0, 0.66
                ] / 1.044;

scatter3(probeLocations(:,1), probeLocations(:,2), probeLocations(:,3), 50, 'r');
plot3(textLocations(1:2,1), textLocations(1:2,2), textLocations(1:2,3), 'r');
plot3(textLocations(3:4,1), textLocations(3:4,2), textLocations(3:4,3), 'r');
text(textLocations(2,1) + 0.1, textLocations(2,2) - 0.1, textLocations(2,3), 'Location A', 'fontName', 'LM Roman 10');
text(textLocations(3,1) + 0.1, textLocations(3,2) + 0.1, textLocations(3,3), 'Location B', 'fontName', 'LM Roman 10');
text(-4.704, 0, 0.66, 'Inlet', 'fontName', 'LM Roman 10', 'rotation', 90, 'horizontalAlignment', 'center');
text(4.704, 0, 0.66, 'Outlet', 'fontName', 'LM Roman 10', 'rotation', 270, 'horizontalAlignment', 'center');

trisurf(windsorBody, 'edgeColor', modelColour, 'faceColor', modelColour);
trisurf(windsorAxles, 'edgeColor', modelColour, 'faceColor', modelColour);
trisurf(windsorWheels, 'edgeColor', modelColour, 'faceColor', modelColour);

% Figure Formatting
view([0, 90]);
axis off;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top', ...
         'dataAspectRatio', [1, 1, 1]);
hold off;

print(fig, '~/MATLAB/Figures/Probe_Locations', '-dpng', '-r300');

clearvars;