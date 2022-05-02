%% Probe Mesh Generator v1.0

clear variables;
close all;
clc;

fig = 0;
figHold = 0;

cellSize = 8e-3; % Base Mesh Cell Size [m]

disp ('=========================');
disp ('Probe Mesh Generator v1.0');
disp ('=========================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit, Combining Functionality of 'surfaceMeshGenerator.m' and 'volumeMeshGenerator.m'


%% Supported OpenFOAM Case Types

% Run_Test
% Windsor_2022


%% Select Mesh Format

disp('Mesh Format');
disp('------------');

disp(' ');

disp('Compatible Mesh Formats:');
disp('    A: Base Surface Mesh');
disp('    B: Wake Volume Mesh');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Data Format [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp (' ');
disp (' ');


%% Select Type Case

disp('Case Type Selection');
disp('--------------------');

caseTypes = {
             'Run_Test'
             'Windsor_SB_wW_Upstream'
             'Windsor_ST_20D_wW_Upstream'
             'Windsor_RSST_16D_U50_wW_Upstream'
            };

valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Case Type', ...
                             'listString', caseTypes);
    
    if ~valid
        disp('WARNING: No Case Type Selected');
    end
    
end

caseType = caseTypes{index};

disp(' ');

disp(['Case: ', caseType]);

disp (' ');
disp (' ');

%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, precision] = selectGeometry(false);

disp (' ');
disp (' ');


%% Delete Previously Generated Meshes

switch format

    case 'A'
        if exist(['Output/Files/surfaceMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm'], 'file')
            delete(['Output/Files/surfaceMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm']);
        end

    case 'B'
        if exist(['Output/Files/volumeMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm'], 'file')
            delete(['Output/Files/volumeMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm']);
        end

end


%% Mesh Generation

disp('Mesh Generation');
disp('----------------');

disp(' ');

disp('***********');
disp('  Running  ');

tic;
switch format

    case 'A'
        % Identify Model Base
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            
            if round(max(geometry.(parts{i}).vertices(:,1)), precision) == xDims(2)
                break
            end
            
            if i == height(parts)
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end
            
        end
    
        geoPoints = geometry.(parts{i}).vertices;
        basePoints = geoPoints(round(geoPoints(:,1), precision) == xDims(2),:);
        
        xDimsBase = xDims(2);
        yDimsBase = round([min(basePoints(:,2)); max(basePoints(:,2))], precision);
        zDimsBase = round([min(basePoints(:,3)); max(basePoints(:,3))], precision);
        
        % Adjust Cell Size to Fit Base
        cellSizeY = (yDimsBase(2) - yDimsBase(1)) / round((yDimsBase(2) - yDimsBase(1)) / cellSize);
        cellSizeZ = (zDimsBase(2) - zDimsBase(1)) / round((zDimsBase(2) - zDimsBase(1)) / cellSize);

        disp(' ');

        % Generate Mesh
        disp('    Generating Initial Surface Mesh...');

        [y, z] = meshgrid(yDimsBase(1):cellSizeY:yDimsBase(2), zDimsBase(1):cellSizeZ:zDimsBase(2));

        % Convert Mesh to Readable Format        
        meshPoints = zeros(height(y(:)),3);
        meshPoints(:,1) = xDimsBase + 1e-6; % Offset to Prevent Base Intersection
        meshPoints(:,(2:3)) = [y(:), z(:)];

        disp(' ');
        
        % Adhere Mesh to Base Boundaries
        disp('    Adhering Mesh to Base Boundaries...');
        
        basePerim = boundary(basePoints(:,2), basePoints(:,3), 1);
        basePerim = basePoints(basePerim,:);
        
        [indexA, indexB] = inpolygon(meshPoints(:,2), meshPoints(:,3), basePerim(:,2), basePerim(:,3));
        meshPoints = meshPoints([indexA, indexB],:);
        
        disp(' ');

        % Write Mesh
        disp('    Writing Mesh Data...');

        fileID = fopen(['Output/Files/surfaceMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm'], 'w');
        formatSpec = '                (%g %g %g)\n';
        
        for i = 1:height(meshPoints)
            fprintf(fileID, formatSpec, meshPoints(i,1), meshPoints(i,2), meshPoints(i,3));
        end
        
        fclose(fileID);
        
        disp(' ');
        
        disp(['    Mesh Written to: ~/MATLAB/Output/Files/surfaceMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm']);
    
    case 'B'
        % Define Mesh Boundaries
        if contains(caseType, 'Run_Test') || (contains(caseType, 'Windsor') && contains(caseType, 'Upstream'))
            xRange = [-1.00625; 2.304]; % [m]
            yRange = [-0.3545; 0.3545]; % [m]
            zRange = [0.001; 0.499]; % [m]
        end
        
        % Ensure Mesh Includes 'xRange(2)'
        xRange(1) = xRange(2) - (floor((xRange(2) - xRange(1)) / cellSize) * cellSize);
        
        % Ensure Mesh Symmetry About y = 0
        yLim = (floor((2 * yRange(2)) / cellSize) * cellSize) / 2;
        yRange = [-yLim; yLim];
        
        % Ensure Mesh Includes 'zRange(2)'
        zRange(2) = 0.339 + (yRange(2) - 0.1945);
        zRange(1) = zRange(2) - (floor(sum(abs(zRange)) / cellSize) * cellSize);

        disp(' ');

        % Generate Initial Mesh
        disp('    Generating Initial Volume Mesh...');

        [x, y, z] = meshgrid(xRange(1):cellSize:xRange(2), ...
                             yRange(1):cellSize:yRange(2), ...
                             zRange(1):cellSize:zRange(2));
        
        % Convert Mesh to Readable Format
        meshPoints = [x(:), y(:), z(:)];
        
        disp(' ');
        
        % Remove Mesh Points Intersecting Geometry
        disp('    Removing Mesh Points Intersecting Geometry...');
        
        parts = fieldnames(geometry);
        
        for i = 1:height(parts)
            geoPoints = unique(geometry.(parts{i}).vertices, 'rows'); % Unique Points in CAD Model
            delaunayTri = delaunayTriangulation(geoPoints);
            
            intersect = find(~isnan(tsearchn(geoPoints, delaunayTri.ConnectivityList, meshPoints))); % Mesh Points Intersecting CAD Model
            meshPoints(intersect,:) = [];
        end

        disp(' ');
        
        % Write Mesh
        disp('    Writing Mesh Data...');
        
        fileID = fopen(['Output/Files/volumeMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm'], 'w');
        formatSpec = '                (%g %g %g)\n';
        
        for i = 1:height(meshPoints)
            fprintf(fileID, formatSpec, meshPoints(i,1), meshPoints(i,2), meshPoints(i,3));
        end
        
        fclose(fileID);

        disp(' ');
        
        disp(['    Mesh Written to: ~/MATLAB/Output/Files/volumeMesh_', caseType, '_', num2str(cellSize * 1e3), 'mm']);

end
executionTime = toc;

disp(' ');

disp(['    Execution Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  Success  ');
disp('***********');


%% Visualisation

% Figure Setup
fig = fig + 1;
set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', 'Mesh_Visualisation');
set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Figure Plotting
parts = fieldnames(geometry);

for i = 1:height(parts)
    patch('faces', geometry.(parts{i}).faces, ...
          'vertices', geometry.(parts{i}).vertices, ...
          'faceColor', ([128, 128, 128] / 255), ...
          'edgeColor', ([128, 128, 128] / 255), ...
          'lineStyle', 'none');
end

scatter3(meshPoints(:,1), meshPoints(:,2), meshPoints(:,3), 10, ([230, 0, 126] / 255), 'filled');

% Figure Formatting
title(' ', 'color', ([254, 254, 254] / 255));
subtitle(' ');
lightangle(30, 30);
lighting gouraud;
axis off;
box off;
view([30, 30]); % 3D
xticks([]);
yticks([]);
zticks([]);
xlabel([]);
ylabel([]);
zlabel([]);
hold off;
