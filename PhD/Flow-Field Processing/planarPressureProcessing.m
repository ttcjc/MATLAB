%% Planar Pressure Processing v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===============================');
disp ('Planar Pressure Processing v2.0');
disp ('===============================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Support ParaView, Probe and Experimental Planar Base Pressure Data (Only)


%% Select Data Format

disp('Data Format');
disp('------------');

disp(' ');

disp('Compatible Data Formats:');
disp('    A: ParaView Planar Data');
disp('    B: Probe Planar Data');
disp('    C: Experimental Planar Data (Varney)');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Data Format [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp (' ');
disp (' ');


%% Initialisation

switch format

    case 'A'
        [caseName, preData, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialisePVdata('p', normalise);
        
        preData = preData.(cell2mat(fieldnames(preData)));

    case 'B'
        [caseName, preData, timePrecision, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(normalise, nProc);    

    case 'C'
        [caseName, preData, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialiseExpData('p', normalise);
        
        preData = preData.(cell2mat(fieldnames(preData)));

end

disp(' ');
disp(' ');


%% Data Formatting

% Specify Map Boundaries
parts = fieldnames(geometry);
for i = 1:height(parts)
    
    if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
        break
    end
    
    if i == height(parts)
        error('Mismatch Between ''xDims'' and Geometry Bounding Box')
    end
    
end

geoPoints = geometry.(parts{i}).vertices;
basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);

mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
mapPerim = basePoints(mapPerim,:);
basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
basePoly = polybuffer(basePoly, -0.0025, 'jointType', 'square');
mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

if ~all(mapPerim(1,:) == mapPerim(end,:))
    mapPerim = vertcat(mapPerim, mapPerim(1,:)); % Close Boundary
end

clear basePoints basePoly;

xLimsData = xDims(2);
yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];

% Adjust Data Origin
switch format

    case {'A', 'B'}
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            
            if normalise
                preData.position(:,1) = preData.position(:,1) + round((1.325 / 1.044), spacePrecision);
                
                if strcmp(preData.planeOrientation, 'YZ')
                    preData.planePosition = preData.planePosition + round((1.325 / 1.044), spacePrecision);
                end
                
            else
                preData.position(:,1) = preData.position(:,1) + 1.325; %#ok<*UNRCH>
                
                if strcmp(preData.planeOrientation, 'YZ')
                    preData.planePosition = preData.planePosition + 1.325;
                end
                
            end
            
        end
        
end

% Shift Data off Base
preData.position(:,1) = preData.position(:,1) + 1e-3;

% Calculate Pressure Coefficient
switch format
    
    case 'A'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            U = 40; % m/s
            rho = 1.269; % kg/m^3
            pRef = 0 * rho; % Pa
            
            preData.CpMean = (preData.pMean - pRef) / (0.5 * rho * U^2);
        end
        
    case 'B'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            U = 40; % m/s
            rho = 1.269; % kg/m^3
            pRef = 0 * rho; % Pa
            
            preData.CpMean = (preData.pMean - pRef) / (0.5 * rho * U^2);
            preData.Cp = cell(height(preData.time),1);
            preData.CpPrime = preData.Cp;
            
            for i = 1:height(preData.time)
                preData.Cp{i} = (preData.p{i} - pRef) / (0.5 * rho * U^2);
                preData.CpPrime{i} = (preData.pPrime{i} - pRef) / (0.5 * rho * U^2);
            end
            
        end
end

% Perform Blockage Correction
switch format
    
    case 'A'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);

            preData.CpMean = (preData.CpMean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
        end
        
    case 'B'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);

            preData.CpMean = (preData.CpMean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
            
            for i = 1:height(preData.time)
                preData.Cp{i} = (preData.Cp{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
                preData.CpPrime{i} = (preData.CpPrime{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
            end
            
        end
end

% Calculate Centre of Pressure
switch format
    
    case {'A', 'C'}
        preData.CoPmean = zeros(1,3);
        preData.CoPmean(1) = preData.position(1,1);
        
        for i = 1:height(preData.position)
            preData.CoPmean(2) = preData.CoPmean(2) + (preData.CpMean(i) * preData.position(i,2));
            preData.CoPmean(3) = preData.CoPmean(3) + (preData.CpMean(i) * preData.position(i,3));
        end
        
        preData.CoPmean(2) = preData.CoPmean(2) / sum(preData.CpMean);
        preData.CoPmean(3) = preData.CoPmean(3) / sum(preData.CpMean);
        
    case 'B'
        preData.CoPmean = zeros(1,3);
        preData.CoPmean(1) = preData.position(1,1);
        
        for i = 1:height(preData.position)
            preData.CoPmean(2) = preData.CoPmean(2) + (preData.CpMean(i) * preData.position(i,2));
            preData.CoPmean(3) = preData.CoPmean(3) + (preData.CpMean(i) * preData.position(i,3));
        end
        
        preData.CoPmean(2) = preData.CoPmean(2) / sum(preData.CpMean);
        preData.CoPmean(3) = preData.CoPmean(3) / sum(preData.CpMean);
        
        preData.CoP = cell(height(preData.time),1);
        
        for i = 1:height(preData.time)
            preData.CoP{i} = zeros(1,3);
            preData.CoP{i}(1) = preData.position(1,1);
            
            for j = 1:height(preData.position)
                preData.CoP{i}(2) = preData.CoP{i}(2) + (preData.Cp{i}(j) * preData.position(j,2));
                preData.CoP{i}(3) = preData.CoP{i}(3) + (preData.Cp{i}(j) * preData.position(j,3));
            end
            
            preData.CoP{i}(2) = preData.CoP{i}(2) / sum(preData.Cp{i});
            preData.CoP{i}(3) = preData.CoP{i}(3) / sum(preData.Cp{i});
        end
        
end

%% Data Presentation

disp('Data Presentation');
disp('------------------');

switch format
    
    case 'B'
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Plot Instantaneous Data? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotInst = false;
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                plotInst = true;
                
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end
            
        end
        
end

% Define Plot Limits
orientation = 'YZ';
        
if contains(caseName, ["Run_Test", "Windsor", "Varney"])
    xLimsPlot = [0.31875; 1.52725];
    yLimsPlot = [-0.2445; 0.2445];
    zLimsPlot = [0; 0.389];
end

if normalise
    xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
    yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
    zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
end

disp(' ');

disp('Presenting Time-Averaged Pressure Data...');

positionData = preData.position;
cMap = viridis(24);
figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
scalarData = preData.CpMean;
figName = [caseName, '_Time_Averaged_Cp'];
CoM = preData.CoPmean;
figSubtitle = ' ';
cLims = [-0.235; -0.023];

fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                        mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                        CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);

disp(' ');

switch format
    
    case 'B'
        
        if plotInst
            disp('Presenting Instantaneous Pressure Data...');
            
            figHold = fig;
            
            for i = 1:height(preData.time)
                
                if i ~= 1
                    clf(fig)
                    fig = figHold;
                end
                
                scalarData = preData.Cp{i};
                figTime = num2str(preData.time(i), ['%.', num2str(timePrecision), 'f']);
                figName = [caseName, '_Instantaneous_Cp_T', erase(figTime, '.')];
                CoM = preData.CoP{i};
                figSubtitle = [num2str(preData.time(i), ['%.', num2str(timePrecision), 'f']), ' \it{s}'];

                fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                        mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                        CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
            end
            
        end
        
end
