%% Lagrangian Volume Field Generator v2.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===========================');
disp ('Volume Field Generator v2.0');
disp ('===========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Follow Recent Lagrangian Processing Structure Changes


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase('global');

disp(' ');
disp(' ');


%% Lagrangian Data Acquisition

disp('LAGRANGIAN DATA ACQUISITION');
disp('---------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Load Saved Particle Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        disp(' ');
        [particleData, particleProps] = lagrangianData(caseFolder, timeDirs);
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/particleData/*.*', ...
                                         'Select Lagrangian Data');

        if contains(filePath, '/MATLAB/particleData/') 
            variableInfo = who('-file', [filePath, fileName]);
            
            if ismember('particleDataGlobal', variableInfo)
                disp(['    Loading: ', fileName]);
                load([filePath, fileName]);
                particleData = particleDataGlobal;
                particleProps = particlePropsGlobal;
                clearvars particleDataGlobal particlePropsGlobal
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

clearvars variableInfo;

disp(' ');
disp(' ');


%% Field Options

disp('FIELD OPTIONS');
disp('-------------');

disp(' ');
disp('Possible Field Locations:');
disp('    A: Downstream (Wake)');
disp('    A: Downstream (Far-Field)');
disp('    C: Upstream (NYI)');

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
    else
        disp('    WARNING: Invalid Entry');
    end

end

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


%% Field Generation

disp('FIELD GENERATION');
disp('-----------------');

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
            xLims = round([0.31875; 1.075] / 1.044, xPre);
            yLims = round([-0.2445; 0.2445] / 1.044, yPre);
            zLims = round([0; 0.389] / 1.044, zPre);
            
        case 'B'
            xLims = round([0.31875; 2.573] / 1.044, xPre);
            yLims = round([-0.3945; 0.3945] / 1.044, yPre);
            zLims = round([0; 0.539] / 1.044, zPre);
            
        case 'C'
            % Add Support for Upstream Spray
    
    end
    
    for i = 1:height(particleData.time)
        particleData.positionCartesian{i,1}(:,1) = round(particleData.positionCartesian{i,1}(:,1) / 1.044, xPre);
        particleData.positionCartesian{i,1}(:,2) = round(particleData.positionCartesian{i,1}(:,2) / 1.044, yPre);
        particleData.positionCartesian{i,1}(:,3) = round(particleData.positionCartesian{i,1}(:,3) / 1.044, zPre);
    end
    
else
    % Add Support for Future Geometries
end

disp(' ');

% Identify Particles of Interest
disp('    Identifying Particles of Interest');

for i = 1:height(particleData.time)
    disp(['        T = ', num2str(str2double(particleData.time{i,1}), '%.4f'), ' s']);
    
    switch method
        
        case {'A', 'B'}
            index = find(particleData.active{i,1} & ...
                         (particleData.d{i,1} * 1e6) >= minD & ...
                         (particleData.d{i,1} * 1e6) <= maxD & ...
                         particleData.positionCartesian{i,1}(:,1) >= xLims(1) & ...
                         particleData.positionCartesian{i,1}(:,1) <= xLims(2) & ...
                         particleData.positionCartesian{i,1}(:,2) >= yLims(1) & ...
                         particleData.positionCartesian{i,1}(:,2) <= yLims(2) & ...
                         particleData.positionCartesian{i,1}(:,3) >= zLims(1) & ...
                         particleData.positionCartesian{i,1}(:,3) <= zLims(2));
                     
        case 'C'
            % Add Support for Upstream Spray
            
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

% Initialise Volume Field
disp('    Generating Volume Field');

cellSize = 5e-3; %l

cellSizeX = (xLims(2) - xLims(1)) / round((xLims(2) - xLims(1)) / cellSize);
cellSizeY = (yLims(2) - yLims(1)) / round((yLims(2) - yLims(1)) / cellSize);
cellSizeZ = (zLims(2) - zLims(1)) / round((zLims(2) - zLims(1)) / cellSize);

[volumeData.y, volumeData.x, volumeData.z] = meshgrid(yLims(1):cellSizeY:yLims(2), ...
                                                      xLims(1):cellSizeX:xLims(2), ...
                                                      zLims(1):cellSizeZ:zLims(2));
                 
for i = 1:height(timeInst)
    
    if isempty(contaminantData.positionCartesian{timeInst(i,1),1})
        disp(['        T = ', num2str(str2double(contaminantData.time{timeInst(i,1),1}), '%.4f'), ...
              ' s - No Particles Present in Volume (Skipping)']);
        continue;
    end
    
    disp(['        T = ', num2str(str2double(contaminantData.time{timeInst(i,1),1}), '%.4f'), ' s']);
    
    volumeData.time{i,1} = contaminantData.time{timeInst(i,1),1};
    
    volumeData.volumeFraction{i,1} = zeros(size(volumeData.x));
    volumeData.nParticle{i,1} = zeros(size(volumeData.x));
    volumeData.d10{i,1} = zeros(size(volumeData.x)); % Arithmetic Mean
    volumeData.d30{i,1} = zeros(size(volumeData.x)); % Volume Mean
    volumeData.d32{i,1} = zeros(size(volumeData.x)); % Sauter Mean
    
    % Assign Contaminant Data to Volume Nodes
    for j = 1:height(contaminantData.positionCartesian{timeInst(i,1),1})
        xDist = abs(contaminantData.positionCartesian{timeInst(i,1),1}(j,1) - volumeData.x(:,1,1));
        yDist = abs(contaminantData.positionCartesian{timeInst(i,1),1}(j,2) - volumeData.y(1,:,1));
        zDist = abs(contaminantData.positionCartesian{timeInst(i,1),1}(j,3) - volumeData.z(1,1,:));

        index = zeros(3,1); % Array Position of Closest Mesh Node
        [~, index(1,1)] = min(xDist);
        [~, index(2,1)] = min(yDist);
        [~, index(3,1)] = min(zDist);
        
        volumeData.volumeFraction{i,1}(index(1,1), index(2,1), index(3,1)) = volumeData.volumeFraction{i,1}(index(1,1), index(2,1), index(3,1)) + ...
                                                                             (contaminantData.nParticle{timeInst(i,1),1}(j,1) * ...
                                                                             ((1 / 12) * tau * contaminantData.d{timeInst(i,1),1}(j,1)^3));
                                                                         
        volumeData.nParticle{i,1}(index(1,1), index(2,1), index(3,1)) = volumeData.nParticle{i,1}(index(1,1), index(2,1), index(3,1)) + ...
                                                                        contaminantData.nParticle{timeInst(i,1),1}(j,1);
            
        volumeData.d10{i,1}(index(1,1), index(2,1), index(3,1)) = volumeData.d10{i,1}(index(1,1), index(2,1), index(3,1)) + ...
                                                                  (contaminantData.nParticle{timeInst(i,1),1}(j,1) * ...
                                                                  contaminantData.d{timeInst(i,1),1}(j,1));
                                                              
        volumeData.d30{i,1}(index(1,1), index(2,1), index(3,1)) = volumeData.volumeFraction{i,1}(index(1,1), index(2,1), index(3,1));
                                                              
        volumeData.d32{i,1}(index(1,1), index(2,1), index(3,1)) = volumeData.d32{i,1}(index(1,1), index(2,1), index(3,1)) + ...
                                                                  (contaminantData.nParticle{timeInst(i,1),1}(j,1) * ...
                                                                  ((1 / 2) * tau * contaminantData.d{timeInst(i,1),1}(j,1)^2));
                                                                  
                                                              

                                                              
    end
    
    volumeData.d10{i,1} = (volumeData.d10{i,1} ./ volumeData.nParticle{i,1}) * 1e6;
    volumeData.d30{i,1} = ((12 * (volumeData.d30{i,1} ./ volumeData.nParticle{i,1}) / tau).^(1 / 3)) * 1e6;
    volumeData.d32{i,1} = (6 * (volumeData.volumeFraction{i,1} ./ volumeData.d32{i,1})) * 1e6;
    volumeData.volumeFraction{i,1} = volumeData.volumeFraction{i,1} / (cellSizeX * cellSizeY * cellSizeZ);
end

% Clear Redundant Data
clearvars contaminantData;

disp(' ');

% Present Contaminant Map
disp('    Presenting Contaminant Map');

figHold = fig;

for i = 1:height(volumeData.time)
    
    if i ~= 1
        clf(fig)
    end
    
    if isempty(volumeData.volumeFraction{i,1})
        continue;
    end
    
    disp(['        T = ', num2str(str2double(volumeData.time{i,1}), '%.4f'), ' s']);
    
    % Plot Volume Field
    switch method
        
        case {'A', 'B'}
            x = volumeData.x;
            y = volumeData.y;
            z = volumeData.z;
            fieldData = volumeData.volumeFraction{i,1};
            fig = figHold;
            isoValue = 5e-6;
            figTime = num2str(str2double(volumeData.time{i,1}), '%.4f');
            figName = ['Spray_Volume_Fraction_', num2str(isoValue), '_T' ...
                       erase(figTime, '.'), '_D', num2str(minD), '_D', num2str(maxD)];
            figTitle = {' ', ' '};
            xLimsPlot = xLims;
            yLimsPlot = yLims;
            zLimsPlot = zLims;
            
            fig = volumeFieldPlots(x, y, z, fieldData, ...
                                   fig, figName, geometry, isoValue, ...
                                   figTitle, xLimsPlot, yLimsPlot, zLimsPlot);
            
        case 'C'
            % Add Support for Upstream Spray
            
    end
    
end

% Format volumeData to Be Human Readable
volumeData.positionGrid = [volumeData.x(:), volumeData.y(:), volumeData.z(:)];

for i = 1:height(volumeData.time)
    volumeData.volumeFraction{i,1} = volumeData.volumeFraction{i,1}(:);
    volumeData.nParticle{i,1} = volumeData.nParticle{i,1}(:);
    volumeData.d10{i,1} = volumeData.d10{i,1}(:);
    volumeData.d30{i,1} = volumeData.d30{i,1}(:);
    volumeData.d32{i,1} = volumeData.d32{i,1}(:);
end

% Clear Redundant Data
volumeData = rmfield(volumeData, {'x', 'y', 'z'});

disp(' ');

executionTime = toc;

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Saving Volume Data

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Volume Data for Future Use? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        namePos = max(strfind(caseFolder, '/')) + 1;

        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/volumeData/', caseFolder(namePos:end)], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/volumeData/', caseFolder(namePos:end)]);
        end

        switch method
            
            case 'A'
                volumeType = 'Wake';

            case 'B'
                volumeType = 'FarField';

            case 'C'
                % Add Support for Upstream Spray

        end

        startInst = erase(num2str(str2double(volumeData.global.time{1,1}), '%.4f'), '.');
        endInst = erase(num2str(str2double(volumeData.global.time{end,1}), '%.4f'), '.');
        
        disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/volumeData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' volumeType, '.mat']);
        save(['/mnt/Processing/Data/Numerical/MATLAB/volumeData/', caseFolder(namePos:end), '/T', startInst, '_T', endInst, '_D', num2str(minD), '_D', num2str(maxD), '_' volumeType, '.mat'], 'mapData', 'minD', 'maxD', '-v7.3', '-noCompression');

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Cleaning

clearvars -except volumeData minD maxD
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