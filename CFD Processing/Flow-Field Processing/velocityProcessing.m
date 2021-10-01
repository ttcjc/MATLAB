%% Velocity Processing v3.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('========================');
disp ('Velocity Processing v3.0');
disp ('========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Improve Versatility
% v3.0 - Moved Plotting to velocityPlots.m


%% Case Initialisation

[caseFolder, data, xDims, yDims, zDims, geometry] = initialisePVdata('U');

disp(' ');
disp(' ');


%% Plotting Options

disp('PLOTTING OPTIONS');
disp('----------------');

disp(' ');
disp('Possible Plotting Regions:');
disp('    A: Near-Wake');
disp('    B: Far-Field (Tight)');
disp('    C: Far-Field (Wide)');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Plotting Region [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        region = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        region = 'B';
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        region = 'C';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');
disp(' ');

% Temporary Implementation Control
if ~contains(caseFolder, 'Windsor')
    error('Case Type Not Yet Supported');
end


%% Velocity Data Processing

disp('PROCESS VELOCITY DATA');
disp('---------------------');
disp(' ');
disp('***********');
disp('  Running  ');

tic;
              
% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

% Flow-Field Limits
if contains(caseFolder, 'Windsor')

    switch region
        
        case 'A'
            xLims = round([0.31875; 1.075] / 1.044, xPre);
            yLims = round([-0.2445; 0.2445] / 1.044, yPre);
            zLims = round([0; 0.489] / 1.044, zPre);
            
        case 'B'
            xLims = round([0.31875; 2.573] / 1.044, xPre);
            yLims = round([-0.2445; 0.2445] / 1.044, yPre);
            zLims = round([0; 0.489] / 1.044, zPre);
        
        case 'C'
            xLims = round([0.31875; 2.573] / 1.044, xPre);
            yLims = round([-0.3945; 0.3945] / 1.044, yPre);
            zLims = round([0; 0.539] / 1.044, zPre);
    
    end
    
end

% Format Plane Data
for i = 1:height(velData.files)
    plane = velData.files{i,1}(1:(end - 4));
    disp(['    Formatting ', plane, ' Data']);
    
    % Store Data
    import = readmatrix([caseFolder, '/', velData.files{i,1}]);
            
    % Shift Data Origin
    if contains(caseFolder, 'Upstream')
        import(:,1) = import(:,1) + 1.325;
    else
        % Add Support for Future Configurations
    end
    
    % Normalise Data
    if contains(caseFolder, 'Windsor')
        velData.(plane).position(:,1) = round(import(:,1) / 1.044, xPre);
        velData.(plane).position(:,2) = round(import(:,2) / 1.044, yPre);
        velData.(plane).position(:,3) = round(import(:,3) / 1.044, zPre);
        velData.(plane).velocity = import(:,4:6) / 40;
    else
        % Add Support for Future Geometries
    end
    
    % Remove Duplicate Entries
    [velData.(plane).position, index, ~] = unique(velData.(plane).position, 'rows', 'stable');
    velData.(plane).velocity = velData.(plane).velocity(index,:);
    
    % Identify Plane Orientation
    if height(unique(velData.(plane).position(:,1))) == 1
        velData.(plane).planeOrientation = 'X';
        velData.(plane).planePosition = unique(velData.(plane).position(:,1));
    elseif height(unique(velData.(plane).position(:,2))) == 1
        velData.(plane).planeOrientation = 'Y';
        velData.(plane).planePosition = unique(velData.(plane).position(:,2));
    elseif height(unique(velData.(plane).position(:,3))) == 1
        velData.(plane).planeOrientation = 'Z';
        velData.(plane).planePosition = unique(velData.(plane).position(:,3));
    else
        error('Invalid Dataset (Unsupported Plane Orientation)');
    end
    
end

disp(' ');

% Present Plane Data
for i = 1:height(velData.files)
    plane = velData.files{i,1}(1:(end - 4));
    disp(['    Presenting ', plane, ' Data']);
    
    planeOrientation = velData.(plane).planeOrientation;
    planePosition = velData.(plane).planePosition;
    positionData = velData.(plane).position;
    velocityData = velData.(plane).velocity;
    
    if contains(planeOrientation, 'X')
        nComponents = 3;
    else
        nComponents = 2;
    end
    
    component = [];
    streamlines = true;
    
    if contains(planeOrientation, 'X')
   
        if contains(caseFolder, 'Windsor')
            modelOutline = [
                            planePosition, -0.1945, 0.004;
                            planePosition, -0.1945, 0.339;
                            planePosition, 0.1945, 0.339;
                            planePosition, 0.1945, 0.004;
                            planePosition, 0.1395, 0.004;
                            planePosition, 0.1395, 0.05;
                            planePosition, -0.1395 0.05;
                            planePosition, -0.1395, 0.004;
                            planePosition, -0.1945, 0.004;
                           ];
        else
            modelOutline = [];
        end
        
    else
        modelOutline = [];
    end
    
    cMap = viridis(24);
    cLims = [0, 1];
    figName = velData.files{i,1}(1:(end - 4));
    figTitle = {' ', ' '};
    
    fig = velocityPlots(planeOrientation, planePosition, ...
                        xLims, yLims, zLims, ...
                        positionData, velocityData, ...
                        nComponents, component, ...
                        fig, figName, cMap, geometry, ...
                        streamlines, modelOutline, figTitle, cLims, ...
                        xLims, yLims, zLims);
end

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Cleaning

clearvars -except data
disp(' ');
