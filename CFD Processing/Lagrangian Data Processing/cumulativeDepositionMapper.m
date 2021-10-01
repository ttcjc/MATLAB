%% Cumulative Lagrangian Deposition Mapper v2.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=================================');
disp ('Cumulative Deposition Mapper v2.0');
disp ('=================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination Only)
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewrite to Improve Versatility and Support Far-Field Extraction Plane


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase('global');

disp(' ');
disp(' ');


%% Mapping Location

disp('MAPPING LOCATION');
disp('----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base-Only');
disp('    B: Full-Body (NYI)');
disp('    C: Far-Field Extraction Plane');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        method = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        method = 'B';
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        method = 'C';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');
disp(' ');

% Temporary Implementation Control
if method == 'B' 
    error('Mapping Location Not Yet Implemented');
end


%% Lagrangian Data Acquisition

disp('LAGRANGIAN DATA ACQUISITION');
disp('---------------------------');

switch method
    
    case {'A', 'B'}
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Load Saved Particle Data? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                disp(' ');
                [particleData, particleProps] = lagrangianDataGlobal(caseFolder, timeDirs);
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/particleData/*.*', ...
                                                 'Select Lagrangian Data');                             

                if contains(filePath, '/MATLAB/particleData/')
                    variableInfo = who('-file', [filePath, fileName]);
                    
                    if ismember('particleDataGlobal', variableInfo)
                        disp(['    Loading: ', fileName]);
                        load([filePath, fileName], 'particleDataGlobal', 'particlePropsGlobal');
                        particleData = particleDataGlobal;
                        particleProps = particlePropsGlobal;
                        clearvars particleDataGlobal particlePropsGlobal;
                        valid = true;
                    else
                        disp('    WARNING: Global Lagrangian Data Not Available in Specified File');
                        clearvars fileName filePath;
                    end
                    
                else
                    disp('    WARNING: Invalid File Selection');
                    clearvars fileName filePath;
                end

            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
    case 'C'
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Load Saved Particle Data? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                disp(' ');
                [particleData, particleProps] = lagrangianDataPlanar(caseFolder, timeDirs);
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/particleData/*.*', ...
                                                 'Select Lagrangian Data');                             
                variableInfo = who('-file', [filePath, fileName]);

                if contains(filePath, '/MATLAB/particleData/')
                    
                    if ismember('particleDataPlanar', variableInfo)
                        disp(['    Loading: ', fileName]);
                        load([filePath, fileName], 'particleDataPlanar', 'particlePropsPlanar');
                        particleData = particleDataPlanar;
                        particleProps = particlePropsPlanar;
                        clearvars particleDataPlanar particlePropsPlanar;
                        valid = true;
                    else
                        disp('    WARNING: Planar Lagrangian Data Not Available in Specified File');
                        clearvars fileName filePath;
                    end
                    
                else
                    disp('    WARNING: Invalid File Selection');
                    clearvars fileName filePath;
                end

            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
end

clearvars variableInfo;

disp(' ');
disp(' ');


%% Mapping Options

disp('MAPPING OPTIONS');
disp('---------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        minD = floor(min(particleData.d{end,1}) * 1e6);
        maxD = ceil(max(particleData.d{end,1}) * 1e6);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        minD = inputD('Minimum');
        maxD = inputD('Maximum');

        if maxD < floor(min(particleData.d{end,1}) * 1e6) || minD > ceil(max(particleData.d{end,1}) * 1e6)
            disp('        WARNING: No Lagrangian Data in Selected Data Range');
        elseif maxD < minD
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');

disp('Data Available for the Following Time Instances:');

for i = 1:height(particleData.time)
    disp(['    ', num2str(particleData.time{i,1}), ' s']);
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Record All Time Instances? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        timeInst = inputTimes(particleData.time);
        
        if isempty(timeInst)
            disp('        WARNING: No Lagrangian Data in Selected Data Range');
        else
            valid = true;
        end
        
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        timeInst = (1:height(particleData.time))';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

disp(' ');
disp(' ');


%% Contaminant Mapping

disp('CONTAMINANT MAPPING');
disp('-------------------');

disp(' ');

disp('***********');
disp('  Running  ');

tic;

disp(' ');

disp('    Initialising');

% Shift Data Origin
if contains(caseFolder, 'Upstream')
    xDims = xDims + 1.325;
    
    for i = 1:height(particleData.time)
        particleData.positionCartesian{i,1}(:,1) = particleData.positionCartesian{i,1}(:,1) + 1.325;
    end
    
else
    % Add Support for Future Configurations
end

% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

if contains(caseFolder, 'Test_Block') || contains(caseFolder, 'Windsor')    
    xDims = round(xDims / 1.044, xPre);
    yDims = round(yDims / 1.044, yPre);
    zDims = round(zDims / 1.044, zPre);
    
    switch method
        
        case 'A'
            xLims = round([-0.61075; 0.53325] / 1.044, xPre);
            yLims = round([-0.2445; 0.2445] / 1.044, yPre);
            zLims = round([0; 0.389] / 1.044, zPre);
            
        case 'B'
            % Add Support for Full-Body Mapping
            
        case 'C'
            planePosition = min(particleData.positionCartesian{1,1}(:,1));
            
            xLims = round([-0.96075; planePosition] / 1.044, xPre);
            yLims = round([-0.5945; 0.5945] / 1.044, yPre);
            zLims = round([0; 0.789] / 1.044, zPre);
            
            planePosition = round(planePosition / 1.044, xPre);
            
    end
    
    for i = 1:height(particleData.time)
        particleData.positionCartesian{i,1}(:,1) = round(particleData.positionCartesian{i,1}(:,1) / 1.044, xPre);
        particleData.positionCartesian{i,1}(:,2) = round(particleData.positionCartesian{i,1}(:,2) / 1.044, yPre);
        particleData.positionCartesian{i,1}(:,3) = round(particleData.positionCartesian{i,1}(:,3) / 1.044, zPre);
    end
    
else
    % Add Support for Future Geometries
end

% Identify Model Boundaries
switch method
    
    case 'A'
        part = fieldnames(geometry);
        for i = 1:height(part)

            if round(max(geometry.(part{i,1}).vertices(:,1)), xPre) == xDims(2)
                break
            end

            if i == height(part)
                disp(' ');
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end

        end

        geoPoints = unique(geometry.(part{i,1}).vertices, 'stable', 'rows');
        index = find(round(geoPoints(:,1), xPre) == xDims(2));

        xDimsBase = xDims(2);
        yDimsBase = round([min(geoPoints(index,2)); max(geoPoints(index,2))], yPre);
        zDimsBase = round([min(geoPoints(index,3)); max(geoPoints(index,3))], zPre);
        
        modelOutline = [];
        
    case 'B'
        % Add Support for Full-Body Mapping
        
    case 'C'
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

            modelOutline(:,2) = round(modelOutline(:,2) / 1.044, yPre);
            modelOutline(:,3) = round(modelOutline(:,3) / 1.044, zPre);
        else
            % Add Support for Future Geometries
            modelOutline = [];
        end

        
end

disp(' ');

% Identify Particles of Interest
disp('    Identifying Particles of Interest');

for i = 1:height(particleData.time)
    disp(['        T = ', num2str(str2double(particleData.time{i,1}), '%.4f'), ' s']);
    
    
    switch method
        
        case 'A'
            index = find(~particleData.active{i,1} & ...
                         (particleData.d{i,1} * 1e6) >= minD & ...
                         (particleData.d{i,1} * 1e6) <= maxD & ...
                         particleData.positionCartesian{i,1}(:,1) == xDimsBase & ...
                         particleData.positionCartesian{i,1}(:,2) >= yDimsBase(1) & ...
                         particleData.positionCartesian{i,1}(:,2) <= yDimsBase(2) & ...
                         particleData.positionCartesian{i,1}(:,3) >= zDimsBase(1) & ...
                         particleData.positionCartesian{i,1}(:,3) <= zDimsBase(2));
                     
        case 'B'
            % Add Support for Full-Body Mapping
                     
        case 'C'
            index = find((particleData.d{i,1} * 1e6) >= minD & ...
                         (particleData.d{i,1} * 1e6) <= maxD & ...
                         particleData.positionCartesian{i,1}(:,2) >= yLims(1) & ...
                         particleData.positionCartesian{i,1}(:,2) <= yLims(2) & ...
                         particleData.positionCartesian{i,1}(:,3) >= zLims(1) & ...
                         particleData.positionCartesian{i,1}(:,3) <= zLims(2));
                     
    end
    
    % Collate Particles of Interest
    contaminantData.time{i,1} = particleData.time{i,1};
    
    for j = 1:height(particleProps)
        prop = particleProps{j,1};
        contaminantData.(prop){i,1} = particleData.(prop){i,1}(index,:);
    end
    
end
    
% Clear Redundant Data
clearvars particleData;

disp(' ');

% Initialise Contaminant Map
disp('    Mapping Contaminants to Grid');

cellSize = 0.01; % l

switch method
    
    case 'A'    
        yDimsBase(1) = yDimsBase(1) + 0.002;
        yDimsBase(2) = yDimsBase(2) - 0.002;
        zDimsBase(1) = zDimsBase(1) + 0.002;
        zDimsBase(2) = zDimsBase(2) - 0.002;
        
        cellSizeX = cellSize;
        cellSizeY = (yDimsBase(2) - yDimsBase(1)) / round((yDimsBase(2) - yDimsBase(1)) / cellSize);
        cellSizeZ = (zDimsBase(2) - zDimsBase(1)) / round((zDimsBase(2) - zDimsBase(1)) / cellSize);
        
        [y, z] = meshgrid(yDimsBase(1):cellSizeY:yDimsBase(2), ...
                          zDimsBase(1):cellSizeZ:zDimsBase(2));
                      
        mapData.global.positionGrid = zeros(height(y(:)),3);
        mapData.global.positionGrid(:,1) = xDimsBase + 0.002;
        mapData.global.positionGrid(:,2:3) = [y(:), z(:)];
        
    case 'B'
        % Add Support for Full-Body Mapping
    
    case 'C'
        cellSizeX = cellSize;
        cellSizeY = (yLims(2) - yLims(1)) / round((yLims(2) - yLims(1)) / cellSize);
        cellSizeZ = (zLims(2) - zLims(1)) / round((zLims(2) - zLims(1)) / cellSize);
        
        [y, z] = meshgrid(yLims(1):cellSizeY:yLims(2), ...
                          zLims(1):cellSizeZ:zLims(2));
                      
        mapData.global.positionGrid = zeros(height(y(:)),3);
        mapData.global.positionGrid(:,1) = planePosition;
        mapData.global.positionGrid(:,2:3) = [y(:), z(:)];
    
end

% Map Contaminants to 'Local' Grid
% - Only Existing Impact Sites Loaded in Grid
for i = 1:height(timeInst)
    
    if isempty(contaminantData.positionCartesian{timeInst(i,1),1})
        disp(['        T = ', num2str(str2double(contaminantData.time{timeInst(i,1),1}), '%.4f'), ...
             ' s - No Particles Recorded on Surface (Skipping)']);
        continue;
    end
    
    disp(['        T = ', num2str(str2double(contaminantData.time{timeInst(i,1),1}), '%.4f'), ' s']);
    
    mapData.local.time{i,1} = contaminantData.time{timeInst(i,1),1};
    
    switch method
        
        case 'A'          
            % Assign Contaminants to Map Nodes
            mapData.tmp.positionCartesian{i,1} = contaminantData.positionCartesian{timeInst(i,1),1};
            mapData.tmp.nParticle{i,1} = contaminantData.nParticle{timeInst(i,1),1};
            mapData.tmp.d{i,1} = contaminantData.d{timeInst(i,1),1};
            
            mapData.tmp.positionCartesian{i,1}(:,1) = xDimsBase + 0.002;
            index = dsearchn(mapData.global.positionGrid, mapData.tmp.positionCartesian{i,1});
            
            for j = 1:height(mapData.tmp.positionCartesian{i,1})
                mapData.tmp.positionCartesian{i,1}(j,:) = mapData.global.positionGrid(index(j,1),:);
            end
            
            mapData.local.occupiedCells{i,1} = unique(mapData.tmp.positionCartesian{i,1}, 'stable', 'rows');
            
        case 'B'
            % Add Support for Full-Body Mapping
            
        case 'C'
            % Assign Contaminants to Map Nodes
            mapData.tmp.positionCartesian{i,1} = vertcat(contaminantData.positionCartesian{(1:timeInst(i,1)),1});
            mapData.tmp.nParticle{i,1} = vertcat(contaminantData.nParticle{(1:timeInst(i,1)),1});
            mapData.tmp.d{i,1} = vertcat(contaminantData.d{(1:timeInst(i,1)),1});
            
            mapData.tmp.positionCartesian{i,1}(:,1) = planePosition;
            index = dsearchn(mapData.global.positionGrid, mapData.tmp.positionCartesian{i,1});
            
            for j = 1:height(mapData.tmp.positionCartesian{i,1})
                mapData.tmp.positionCartesian{i,1}(j,:) = mapData.global.positionGrid(index(j,1),:);
            end
            
            mapData.local.occupiedCells{i,1} = unique(mapData.tmp.positionCartesian{i,1}, 'stable', 'rows');
    
    end
    
    % Calculate Contaminant Variables in Map Nodes
    mapData.tmp.volume{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1);
    
    mapData.local.nParticle{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1);
    mapData.local.mass{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1);
    mapData.local.d10{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1); % Arithmetic Mean
    mapData.local.d30{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1); % Volume Mean
    mapData.local.d32{i,1} = zeros(height(mapData.local.occupiedCells{i,1}),1); % Sauter Mean
    
    for j = 1:height(mapData.local.occupiedCells{i,1})
        index = find(ismember(mapData.tmp.positionCartesian{i,1}, mapData.local.occupiedCells{i,1}(j,:), 'rows'));
        
        for k = 1:height(index)
            mapData.tmp.volume{i,1}(j,1) = mapData.tmp.volume{i,1}(j,1) + ...
                                           (mapData.tmp.nParticle{i,1}(index(k,1),1) * ...
                                           ((1 / 12) * tau * mapData.tmp.d{i,1}(index(k,1),1)^3));
            
            mapData.local.nParticle{i,1}(j,1) = mapData.local.nParticle{i,1}(j,1) + ...
                                                mapData.tmp.nParticle{i,1}(index(k,1),1);
                                            
            mapData.local.d10{i,1}(j,1) = mapData.local.d10{i,1}(j,1) + ...
                                          (mapData.tmp.nParticle{i,1}(index(k,1),1) * ...
                                          mapData.tmp.d{i,1}(index(k,1),1));
                                      
            mapData.local.d32{i,1}(j,1) = mapData.local.d32{i,1}(j,1) + ...
                                          (mapData.tmp.nParticle{i,1}(index(k,1),1) * ...
                                          ((1 / 2) * tau * mapData.tmp.d{i,1}(index(k,1),1)^2));
        end
        
    end
    
    mapData.local.mass{i,1} = 1000 * mapData.tmp.volume{i,1};
    mapData.local.massNorm{i,1} = mapData.local.mass{i,1} / 2.933392467046320e-07; % Square-Back Base Peak
    mapData.local.d10{i,1} = (mapData.local.d10{i,1} ./ mapData.local.nParticle{i,1}) * 1e6;
    mapData.local.d30{i,1} = ((12 * (mapData.tmp.volume{i,1} ./ mapData.local.nParticle{i,1}) / tau).^(1 / 3)) * 1e6;
    mapData.local.d32{i,1} = (6 * (mapData.tmp.volume{i,1} ./ mapData.local.d32{i,1})) * 1e6;    
end

% Clear Redundant Data
clearvars contaminantData;
mapData = rmfield(mapData, 'tmp');

i = 1;
while i <= height(mapData.local.occupiedCells)
    
    if isempty(mapData.local.occupiedCells{i,1})
        mapData.local.time(i,:) = [];
        mapData.local.occupiedCells(i,:) = [];
        mapData.local.nParticle(i,:) = [];
        mapData.local.mass(i,:) = [];
        mapData.local.massNorm(i,:) = [];
        mapData.local.d10(i,:) = [];
        mapData.local.d30(i,:) = [];
        mapData.local.d32(i,:) = [];
    else
        i = i + 1;
    end
    
end

disp(' ');

% Restructure Grid to Span All Time Instances (Global Grid)
disp('    Restructuring Contaminant Grid');

for i = 1:height(mapData.local.time)
    
    if isempty(mapData.local.occupiedCells{i,1})
        continue;
    end
    
    disp(['        T = ', num2str(str2double(mapData.local.time{i,1}), '%.4f'), ' s']);

    mapData.global.time{i,1} = mapData.local.time{i,1};
    mapData.global.nParticle{i,1} = zeros(height(mapData.global.positionGrid),1);
    mapData.global.mass{i,1} = zeros(height(mapData.global.positionGrid),1);
    mapData.global.massNorm{i,1} = zeros(height(mapData.global.positionGrid),1);
    mapData.global.d10{i,1} = zeros(height(mapData.global.positionGrid),1);
    mapData.global.d30{i,1} = zeros(height(mapData.global.positionGrid),1);
    mapData.global.d32{i,1} = zeros(height(mapData.global.positionGrid),1);
    
    for j = 1:height(mapData.local.occupiedCells{i,1})
        index = find(ismember(mapData.global.positionGrid, mapData.local.occupiedCells{i,1}(j,:), 'rows'));
        
        mapData.global.nParticle{i,1}(index,:) = mapData.local.nParticle{i,1}(j,1);
        mapData.global.mass{i,1}(index,:) = mapData.local.mass{i,1}(j,1);
        mapData.global.massNorm{i,1}(index,:) = mapData.local.massNorm{i,1}(j,1);
        mapData.global.d10{i,1}(index,:) = mapData.local.d10{i,1}(j,1);
        mapData.global.d30{i,1}(index,:) = mapData.local.d30{i,1}(j,1);
        mapData.global.d32{i,1}(index,:) = mapData.local.d32{i,1}(j,1);
    end
    
end

disp(' ');

% Calculate Centre of Mass
disp('    Calculating Centre of Mass');

for i = 1:height(mapData.global.time)
    
    if isempty(mapData.global.mass{i,1})
        continue;
    end
    
    disp(['        T = ', num2str(str2double(mapData.global.time{i,1}), '%.4f'), ' s']);
    
    mapData.global.CoM{i,1} = zeros(1,3);
    
    switch method
        
        case 'A'
            mapData.global.CoM{i,1}(1,1) = xDimsBase + 0.002;
            mapData.global.CoM{i,1}(1,2) = sum(mapData.global.mass{i,1}(:,1) .* mapData.global.positionGrid(:,2)) / ...
                                           sum(mapData.global.mass{i,1}(:,1));
            mapData.global.CoM{i,1}(1,3) = sum(mapData.global.mass{i,1}(:,1) .* mapData.global.positionGrid(:,3)) / ...
                                           sum(mapData.global.mass{i,1}(:,1));
                                       
        case 'B'
            % Add Support for Full-Body Mapping
        
        case 'C'
            mapData.global.CoM{i,1}(1,1) = planePosition;
            mapData.global.CoM{i,1}(1,2) = sum(mapData.global.mass{i,1}(:,1) .* mapData.global.positionGrid(:,2)) / ...
                                           sum(mapData.global.mass{i,1}(:,1));
            mapData.global.CoM{i,1}(1,3) = sum(mapData.global.mass{i,1}(:,1) .* mapData.global.positionGrid(:,3)) / ...
                                           sum(mapData.global.mass{i,1}(:,1));
                                       
    end
    
end

disp(' ');

% Present Contaminant Map
disp('    Presenting Contaminant Map');

figHold = fig;

for i = 1:height(mapData.global.time)
    
    if i ~= 1
        clf(fig);
    end
    
    if isempty(mapData.global.massNorm{i,1})
        continue;
    end
    
    disp(['        T = ', num2str(str2double(mapData.global.time{i,1}), '%.4f'), ' s']);
    
    % Plot Contaminant Map
    switch method
        
        case 'A'
            xLimsData = xDimsBase + 0.002;
            yLimsData = yDimsBase;
            zLimsData = zDimsBase;
            positionData = mapData.global.positionGrid;
            contaminantData = mapData.global.massNorm{i,1};
            fig = figHold;
            figTime = num2str(str2double(mapData.global.time{i,1}), '%.4f');
            figName = ['Cumulative_Base_Contamination_Map_T', ...
                       erase(figTime, '.'), '_D', num2str(minD), '_D', num2str(maxD)];
            cMap = viridis(24);
            CoM = mapData.global.CoM{i,1};
            
            figTitle = {' ', ' '};
%             figTitle = {'Time (\it{s})', figTime};
            
            cLims = [0, 1]; % Windsor Full Duration
            xLimsPlot = xLims;
            yLimsPlot = yLims;
            zLimsPlot = zLims;
            
            fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                   positionData, contaminantData, ...
                                   fig, figName, cMap, geometry, ...
                                   CoM, modelOutline, figTitle, cLims, ...
                                   xLimsPlot, yLimsPlot, zLimsPlot);
                       
        case 'B'
            % Add Support for Full-Body Mapping
            
        case 'C'
            xLimsData = planePosition;
            yLimsData = yLims;
            zLimsData = zLims;
            positionData = mapData.global.positionGrid;
            contaminantData = mapData.global.d32{i,1};
            fig = figHold;
            figTime = num2str(str2double(mapData.global.time{i,1}), '%.4f');
            figName = ['Cumulative_Planar_Contamination_Map_T', ...
                       erase(figTime, '.'), '_D', num2str(minD), '_D', num2str(maxD)];
            cMap = viridis(24);
            CoM = mapData.global.CoM{i,1};
            
            figTitle = {' ', ' '};
%             figTitle = {'Time (\it{s})', figTime};
            
%             cLims = [0, 12]; % Windsor Full Duration
            cLims = [0, 120];
            xLimsPlot = xLims;
            yLimsPlot = yLims;
            zLimsPlot = zLims;
            
            fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                positionData, contaminantData, ...
                                fig, figName, cMap, geometry, ...
                                CoM, modelOutline, figTitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot);
                       
    end
                       
end

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Saving Map Data

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Map Data for Future Use? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        namePos = max(strfind(caseFolder, '/')) + 1;

        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/mapDataCumulative/', caseFolder(namePos:end)], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/mapDataCumulative/', caseFolder(namePos:end)]);
        end

        switch method

            case 'A'
                mapType = 'Base';

            case 'B'
                % Add Support for Full-Body Mapping

            case 'C'
                mapType = 'Planar';

        end

        startInst = erase(num2str(str2double(mapData.global.time{1,1}), '%.4f'), '.');
        endInst = erase(num2str(str2double(mapData.global.time{end,1}), '%.4f'), '.');
        
        disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/mapDataCumulative/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' mapType, '.mat']);
        save(['/mnt/Processing/Data/Numerical/MATLAB/mapDataCumulative/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' mapType, '.mat'], 'mapData', 'minD', 'maxD', '-v7.3', '-noCompression');

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Cleaning

clearvars -except mapData minD maxD
disp(' ');


%% Local Functions

function D = inputD(type)

    valid = false;
    while ~valid
        D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));

        if isnan(D) || length(D) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end


function T = inputTimes(time)

    valid = false;
    while ~valid
        T = str2num(input('    List Desired Time Instances (Row Vector Form): ', 's')); %#ok<ST2NM>

        if any(isnan(T)) || ~isrow(T)
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

    T = find(ismember(time, strsplit(num2str(T))));

end
