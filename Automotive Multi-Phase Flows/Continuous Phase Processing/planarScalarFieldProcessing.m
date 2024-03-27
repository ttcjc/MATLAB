%% Planar Scalar Field Processing v1.1
% ----
% Load, Process and Present Numerical Planar Scalar Fields


%% Preamble

run preamble;

normDims = true; % Normalise Spatial Dimensions

figSave = false; % Save .fig File(s)

fieldList = {
             'p'
             'CpT'
             'kResolved'
            };

disp('===================================');
disp('Planar Scalar Field Processing v1.1');
disp('===================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Update To Support Total Pressure Coefficient


%% Select Variable of Interest

disp('Select Variable to Present...');

valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Scalar Field', ...
                             'listString', fieldList);

    if ~valid
        disp('    WARNING: No Mapping Variable Selected');
    end
    
end
clear valid;

field = fieldList{index}; clear fieldVars;

disp(['    Field of Interest: ', field]);

disp(' ');
disp(' ');


%% Acquire Field Data

[campaignID, caseID, fieldData] = initialisePVdata(saveLoc, field);

planes = fieldnames(fieldData);

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');


%% Process Field Data

disp('Field Data Processing');
disp('----------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

%%%%

disp(' ');

% Adjust Data Origin
disp('    Adjusting Data Origin...');

if strcmp(campaignID, 'Windsor_Upstream_2023')

    for i = 1:height(planes)
        fieldData.(planes{i}).positionGrid(:,1) = fieldData.(planes{i}).positionGrid(:,1) + 1.325;
    end
    clear i;

end

for i = 1:height(planes)
    [fieldData.(planes{i}).positionGrid, index] = unique(fieldData.(planes{i}).positionGrid, 'rows', 'stable');

    fieldData.(planes{i}).(field).mean = fieldData.(planes{i}).(field).mean(index);
end
clear i;

disp(' ');

% Map Raw Data Onto Uniform Grid
disp('    Mapping Raw Data Onto Uniform Grid...');

% Set Target Spatial Resolution
if strcmp(campaignID, 'Windsor_fullScale')
    cellSize.target = 4e-3;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    cellSize.target = 1e-3;
else
    cellSize.target = 1e-3;
end

for i = 1:height(planes)

    if contains(planes{i}, '_X_')
        orientation = 'YZ';
    elseif contains(planes{i}, '_Y_')
        orientation = 'XZ';
    else
        orientation = 'XY';
    end

    switch orientation

        case 'YZ'
            xLimsData = fieldData.(planes{i}).positionGrid(1,1);
            yLimsData = [min(fieldData.(planes{i}).positionGrid(:,2)); ...
                         max(fieldData.(planes{i}).positionGrid(:,2))];
            zLimsData = [min(fieldData.(planes{i}).positionGrid(:,3)); ...
                         max(fieldData.(planes{i}).positionGrid(:,3))];

            % Adjust Uniform Cell Size to Fit Region of Interest
            nPy = (diff(yLimsData) / cellSize.target) + 1;
            nPz = (diff(zLimsData) / cellSize.target) + 1;

            sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
            sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

            cellSize.(planes{i}).x = cellSize.target;
            cellSize.(planes{i}).y = sizeY(1); clear sizeY;
            cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
            cellSize.(planes{i}).area = cellSize.(planes{i}).y * cellSize.(planes{i}).z;

            yOrig = fieldData.(planes{i}).positionGrid(:,2);
            zOrig = fieldData.(planes{i}).positionGrid(:,3);

            % Generate Grid
            [y, z] = ndgrid(linspace(yLimsData(1), yLimsData(2), nPy), ...
                            linspace(zLimsData(1), zLimsData(2), nPz));

            fieldData.(planes{i}).positionGrid = zeros([height(y(:)),3]);
            fieldData.(planes{i}).positionGrid(:,1) = xLimsData;
            fieldData.(planes{i}).positionGrid(:,[2,3]) = [y(:), z(:)];
            
            % Perform Interpolation
            fieldInterp = scatteredInterpolant(yOrig, zOrig, fieldData.(planes{i}).(field).mean, ...
                                               'linear', 'none');

            fieldData.(planes{i}).(field).mean = fieldInterp(fieldData.(planes{i}).positionGrid(:,2), ...
                                                             fieldData.(planes{i}).positionGrid(:,3));

        case 'XZ'
            xLimsData = [min(fieldData.(planes{i}).positionGrid(:,1)); ...
                         max(fieldData.(planes{i}).positionGrid(:,1))];
            yLimsData = fieldData.(planes{i}).positionGrid(1,2);
            zLimsData = [min(fieldData.(planes{i}).positionGrid(:,3)); ...
                         max(fieldData.(planes{i}).positionGrid(:,3))];

            % Adjust Uniform Cell Size to Fit Region of Interest
            nPx = (diff(xLimsData) / cellSize.target) + 1;
            nPz = (diff(zLimsData) / cellSize.target) + 1;

            sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
            sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

            cellSize.(planes{i}).x = sizeX(1); clear sizeX;
            cellSize.(planes{i}).y = cellSize.target;
            cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
            cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).z;

            xOrig = fieldData.(planes{i}).positionGrid(:,1);
            zOrig = fieldData.(planes{i}).positionGrid(:,3);

            % Generate Grid
            [x, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                            linspace(zLimsData(1), zLimsData(2), nPz));

            fieldData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
            fieldData.(planes{i}).positionGrid(:,2) = yLimsData;
            fieldData.(planes{i}).positionGrid(:,[1,3]) = [x(:), z(:)];
            
            % Perform Interpolation
            fieldInterp = scatteredInterpolant(xOrig, zOrig, fieldData.(planes{i}).(field).mean, ...
                                               'linear', 'none');

            fieldData.(planes{i}).(field).mean = fieldInterp(fieldData.(planes{i}).positionGrid(:,1), ...
                                                             fieldData.(planes{i}).positionGrid(:,3));

        case 'XY'
            xLimsData = [min(fieldData.(planes{i}).positionGrid(:,1)); ...
                         max(fieldData.(planes{i}).positionGrid(:,1))];
            yLimsData = [min(fieldData.(planes{i}).positionGrid(:,2)); ...
                         max(fieldData.(planes{i}).positionGrid(:,2))];
            zLimsData = fieldData.(planes{i}).positionGrid(1,3);

            % Adjust Uniform Cell Size to Fit Region of Interest
            nPx = (diff(xLimsData) / cellSize.target) + 1;
            nPy = (diff(yLimsData) / cellSize.target) + 1;

            sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
            sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));

            cellSize.(planes{i}).x = sizeX(1); clear sizeX;
            cellSize.(planes{i}).y = sizeY(1); clear sizeY;
            cellSize.(planes{i}).z = cellSize.target;
            cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).y;

            xOrig = fieldData.(planes{i}).positionGrid(:,1);
            yOrig = fieldData.(planes{i}).positionGrid(:,2);

            % Generate Grid
            [x, y] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                            linspace(yLimsData(1), yLimsData(2), nPy));

            fieldData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
            fieldData.(planes{i}).positionGrid(:,3) = zLimsData;
            fieldData.(planes{i}).positionGrid(:,[1,2]) = [x(:), y(:)];
            
            % Perform Interpolation
            fieldInterp = scatteredInterpolant(xOrig, yOrig, fieldData.(planes{i}).(field).mean, ...
                                               'linear', 'none');

            fieldData.(planes{i}).(field).mean = fieldInterp(fieldData.(planes{i}).positionGrid(:,1), ...
                                                             fieldData.(planes{i}).positionGrid(:,2));

    end
    clear yOrig zOrig y z fieldInterp;

end
clear i;

% Normalise Coordinate System
if normDims
    disp(' ');
    
    disp('    Normalising Spatial Dimensions...');
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
    end
    clear i parts;
    
    xDims = xDims / normLength;
    yDims = yDims / normLength;
    zDims = zDims / normLength;
    
    cellSize.target = cellSize.target / normLength;
    
    for i = 1:height(planes)
        cellSize.(planes{i}).x = cellSize.(planes{i}).x / normLength;
        cellSize.(planes{i}).y = cellSize.(planes{i}).y / normLength;
        cellSize.(planes{i}).z = cellSize.(planes{i}).z / normLength;
        cellSize.(planes{i}).area = cellSize.(planes{i}).area / (normLength^2);

        fieldData.(planes{i}).positionGrid = fieldData.(planes{i}).positionGrid / normLength;
    end
    
end

%%%%

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Map(s)? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMean = false;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMean = true;

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotMean
    
    % Select Plane(s) of Interest
        valid = false;
        while ~valid
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'multiple', ...
                                     'name', 'Select Variable(s) to Plot', ...
                                     'listString', planes);
        
            if ~valid
                disp(    'WARNING: No Planes Selected');
            end
        end
        clear valid;
        
        plotPlanes = planes(index);
end

disp(' ');
disp(' ');

%% Present Pressure Data

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean
    
    if strcmp(campaignID, 'Windsor_fullScale')
        spatialRes = 2e-3;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        spatialRes = 0.5e-3;
    else
        spatialRes = 0.5e-3;
    end
    
    if normDims
        spatialRes = spatialRes / normLength;
    end
    
    if strcmp(field, 'CpT')
        cMap = plasma(32);
        cLims = [-0.8; 1]; % cLims = [-1.2; 1];
    elseif strcmp(field, 'kResolved')
        cMap = cool2warm(32);
        cLims = [0.6; 1];
    else
        cMap = viridis(32);
        cLims = 'auto';
    end
        
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    contourlines = [];
    refPoint = [];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    
    
    for i = 1:height(plotPlanes)
        disp(['Presenting ', plotPlanes{i}, ' Field Data...']);
        
        if contains(plotPlanes{i}, '_X_')
            orientation = 'YZ';
        elseif contains(plotPlanes{i}, '_Y_')
            orientation = 'XZ';
        else
            orientation = 'XY';
        end
        
        positionData = fieldData.(plotPlanes{i}).positionGrid;
        scalarData = fieldData.(plotPlanes{i}).(field).mean;
        
        switch orientation
            
            case 'YZ'
                xLimsPlot = [0.3; 4.6257662];
                yLimsPlot = [-0.5; 0.5];
                zLimsPlot = [0; 0.5];
                
            case {'XZ', 'XY'}
                xLimsPlot = [0.3; 1.2];
                yLimsPlot = [-0.3; 0.3];
                zLimsPlot = [0; 0.5];
                
        end

        if ~normDims
            xLimsPlot = xLimsPlot * normLength;
            yLimsPlot = yLimsPlot * normLength;
            zLimsPlot = zLimsPlot * normLength;
        end
        
        switch orientation

            case 'YZ'
                xLimsData = fieldData.(plotPlanes{i}).positionGrid(1,1);
                yLimsData = yLimsPlot;
                zLimsData = zLimsPlot;

            case 'XZ'
                xLimsData = xLimsPlot;
                yLimsData = fieldData.(plotPlanes{i}).positionGrid(1,2);
                zLimsData = zLimsPlot;

            case 'XY'
                xLimsData = xLimsPlot;
                yLimsData = yLimsPlot;
                zLimsData = fieldData.(plotPlanes{i}).positionGrid(1,3);

        end
        
        figName = [plotPlanes{i}, '_', caseID];
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    
    disp(' ');
end
    
if ~plotMean
    disp('Skipping Map Presentation...');
    
    disp(' ');
end


