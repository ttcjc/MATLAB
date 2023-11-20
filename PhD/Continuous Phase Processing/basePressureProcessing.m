%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Pressure Processing v3.0

figSave = false; % Save .fig File(s)

normDims = true; % Normalise Spatial Dimensions

disp('===============================');
disp('Planar Pressure Processing v3.0');
disp('===============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Support ParaView, Probe and Experimental Base Pressure Data
% v2.1 - Minor Update to Support Changes to 'plotPlanarScalarField'
% v3.0 - Update To Improve Consistency of Structures Across Repository


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
clear valid;

disp(' ');
disp(' ');


%% Acquire Pressure Data

switch format

    case 'A'
        [campaignID, caseID, pData] = initialisePVdata(saveLocation, 'p');
        
        pData = pData.(cell2mat(fieldnames(pData)));
        
    case 'B'
        error('NYI');
        
    case 'C'
        [campaignID, caseID, pData] = initialiseExpFlowData(saveLocation, 'p');
        
        pData = pData.(cell2mat(fieldnames(pData)));
        pData.p = rmfield(pData.p, 'RMS');

end
    
%     case 'B'
%         [caseID, dataID, pData, sampleInterval, timePrecision, geometry, ...
%          xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(saveLocation, normalise, nProc);    

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(normDims);

disp(' ');
disp(' ');

%% Generate Base Pressure Distributions

disp('Base Pressure Distribution Generation');
disp('--------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

%%%%

disp(' ');

% Adjust Data Origin
disp('    Adjusting Data Origin...');

switch format

    case {'A', 'B'}
        
        if ~strcmp(campaignID, 'Windsor_fullScale')
            pData.positionGrid(:,1) = pData.positionGrid(:,1) + 1.325;
        end
        
end

disp(' ');

% Normalise Spatial Dimensions
disp('    Normalising Spatial Dimensions...');

if normDims
    pData.positionGrid = round((pData.positionGrid / normLength), spacePrecision);
end

disp(' ');

% Identify Base Boundaries
disp('    Identifying Base Boundaries...');

parts = fieldnames(geometry);
for i = 1:height(parts)
    
    if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
        parts = parts{i};
        break;
    end
    
    if i == height(parts)
        error('Mismatch Between ''xDims'' and Geometry Bounding Box')
    end
    
end
clear i;

geoPoints = geometry.(parts).vertices;
basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);

mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
mapPerim = basePoints(mapPerim,:);
basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);

if normDims
    
    if strcmp(campaignID, 'Varney')
        basePoly = polybuffer(basePoly, -16e-3, 'jointType', 'square');
    else
        basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
    end
    
else
    
    if strcmp(campaignID, 'Windsor_fullScale') || strcmp(campaignID, 'Varney')
        basePoly = polybuffer(basePoly, -16e-3, 'jointType', 'square');
    else
        basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
    end
    
end

mapPerim = ones([height(basePoly.Vertices),3]) * mapPerim(1,1);
mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

if ~all(mapPerim(1,:) == mapPerim(end,:))
    mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
end

xLimsData = xDims(2);
yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];

disp(' ');

% Map Raw Data Onto Uniform Grid
disp('    Mapping Raw Data Onto Uniform Grid...');

if normDims
    cellSize.target = 1e-3;
else

    if strcmp(campaignID, 'Windsor_fullScale')
        cellSize.target = 4e-3;
    else
        cellSize.target = 1e-3;
    end

end

cellSize.x = cellSize.target;
cellSize.y = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
cellSize.z = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize.target));

cellSize.area = cellSize.y * cellSize.z;

yOrig = pData.positionGrid(:,2);
zOrig = pData.positionGrid(:,3);

[y, z] = ndgrid(yLimsData(1):cellSize.y:yLimsData(2), zLimsData(1):cellSize.z:zLimsData(2));

pData.positionGrid = zeros([height(y(:)),3]);
pData.positionGrid(:,1) = xLimsData;
pData.positionGrid(:,[2,3]) = [y(:), z(:)];

switch format
    
    case 'A'
        pInterp = scatteredInterpolant(yOrig, zOrig, pData.p.mean, ...
                                      'linear', 'none');
        
        pData.p.mean = pInterp(pData.positionGrid(:,2), pData.positionGrid(:,3));
        
    case 'C'
        pInterp = scatteredInterpolant(yOrig, zOrig, pData.Cp.mean, ...
                                      'linear', 'none');
        
        pData.Cp.mean = pInterp(pData.positionGrid(:,2), pData.positionGrid(:,3));
        
        pData.Cp.mean = reshape(pData.Cp.mean, [height(unique(pData.positionGrid(:,2))), ...
                                                height(unique(pData.positionGrid(:,3)))]);
        pData.Cp.mean = smooth2(pData.Cp.mean, 9);
        pData.Cp.mean = pData.Cp.mean(:);
        
end
clear yOrig zOrig y z pInterp;

disp(' ');

% Calculate Pressure Coefficient
disp('    Calculating Pressure Coefficient...');

switch format
    
    case 'A'
        
        if strcmp(campaignID, 'Windsor_fullScale')
            U = 22.22; % m/s
            rho = 1.269; % kg/m^3
            pRef = 0 * rho; % Pa
        else
            U = 40; % m/s
            rho = 1.269; % kg/m^3
            pRef = 19.524 * rho; % Pa
        end
        
        pData.p.mean = pData.p.mean * rho;
        pData.Cp.mean = (pData.p.mean - pRef) / (0.5 * rho * U^2);
        
end
        
%     case 'B'
%         
%         if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
%             U = 40; % m/s
%             rho = 1.269; % kg/m^3
%             pRef = 0 * rho; % Pa
%             
%             pData.CpMean = ((pData.pMean * rho) - pRef) / (0.5 * rho * U^2);
%             pData.Cp = cell(height(pData.time),1);
%             pData.CpPrime = pData.Cp;
%             
%             for i = 1:height(pData.time)
%                 pData.Cp{i} = ((pData.p{i} * rho) - pRef) / (0.5 * rho * U^2);
%                 pData.CpPrime{i} = ((pData.pPrime{i} * rho) - pRef) / (0.5 * rho * U^2);
%             end
%             
%         end

disp(' ');

% Perform Blockage Correction
disp('    Performing Blockage Correction...');

switch format
    
    case 'A'
        
        if ~strcmp(campaignID, 'Windsor_fullScale')
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
            
            pData.p.mean = (pData.p.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
            pData.Cp.mean = (pData.Cp.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
        end
        
end
        
%     case 'B'
%         
%         if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
%             Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
%             At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
% 
%             pData.CpMean = (pData.CpMean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
%             
%             for i = 1:height(pData.time)
%                 pData.Cp{i} = (pData.Cp{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
%                 pData.CpPrime{i} = (pData.CpPrime{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
%             end
%             
%         end

disp(' ');

% Calculate Centre of Pressure
disp('    Calculating Centre of Pressure...');

switch format
    
    case {'A', 'C'}
        pData.CoP.mean = zeros([1,3]);
        
        pData.CoP.mean(1) = pData.positionGrid(1,1);
        pData.CoP.mean(2) = sum(pData.Cp.mean .* pData.positionGrid(:,2)) / sum(pData.Cp.mean);
        pData.CoP.mean(3) = sum(pData.Cp.mean .* pData.positionGrid(:,3)) / sum(pData.Cp.mean);
        
end
        
%     case 'B'
%         pData.CoPmean = zeros([1,3]);
%         
%         pData.CoPmean(1) = pData.position(1,1);
%         pData.CoPmean(2) = sum(pData.CpMean .* pData.position(:,2)) / sum(pData.CpMean);
%         pData.CoPmean(3) = sum(pData.CpMean .* pData.position(:,3)) / sum(pData.CpMean);
%         
%         pData.CoP = cell(height(pData.time),1);
%         
%         for i = 1:height(pData.time)
%             pData.CoP{i} = zeros([1,3]);
%         
%             pData.CoP{i}(1) = pData.position(1,1);
%             pData.CoP{i}(2) = sum(pData.Cp{i} .* pData.position(:,2)) / sum(pData.Cp{i});
%             pData.CoP{i}(3) = sum(pData.Cp{i} .* pData.position(:,3)) / sum(pData.Cp{i});
%         end

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
    
    selection = input('Plot Time-Averaged Map? [y/n]: ', 's');

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

switch format
    
    case 'B'
        valid = false;
        while ~valid
            disp(' ');

            selection = input('Plot RMS Map? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotRMS = false;
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                plotRMS = true;
                
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;

        valid = false;
        while ~valid
            disp(' ');

            selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotInst = false;
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                plotInst = true;

                startFrame = inputFrames(nTimes, 'Start');

                if startFrame == -1
                    continue;
                end

                endFrame = inputFrames(nTimes, 'End');

                if endFrame == -1
                    continue;
                elseif endFrame < startFrame
                    disp('        WARNING: Invalid Time Format (''endFrame'' Precedes ''startFrame'')');

                    continue;
                end

                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;
        
    otherwise
        plotRMS = false;
        plotInst = false;
        
end

disp(' ');
disp(' ');


%% Present Pressure Data

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    
    if normDims
        spatialRes = 0.5e-3;
    else
        
        if strcmp(campaignID, 'Windsor_fullScale')
            spatialRes = 2e-3;
        else
            spatialRes = 0.5e-3;
        end
        
    end
    
    xLimsData = xLimsData + spatialRes;
    positionData = pData.positionGrid; positionData(:,1) = xLimsData;
    nPlanes = 1;
    planeNo = 1;
    cMap = viridis(32);
    contourlines = [];
    
    xLimsPlot = [0.3; 4.6257662];
    yLimsPlot = [-0.25; 0.25];
    zLimsPlot = [0; 0.4];
        
    if ~normDims

        if strcmp(campaignID, 'Windsor_fullScale')
            xLimsPlot = round((xLimsPlot * 4.176), spacePrecision);
            yLimsPlot = round((yLimsPlot * 4.176), spacePrecision);
            zLimsPlot = round((zLimsPlot * 4.176), spacePrecision);
        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            xLimsPlot = round((xLimsPlot * 1.044), spacePrecision);
            yLimsPlot = round((yLimsPlot * 1.044), spacePrecision);
            zLimsPlot = round((zLimsPlot * 1.044), spacePrecision);
        end

    end

    if plotMean
        disp('    Presenting Time-Averaged Pressure Data...');
        
        scalarData = pData.Cp.mean;
        refPoint = [];
        figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
        figName = ['Average_Cp_', caseID];

        if strcmp(campaignID, 'Windsor_fullScale')
            cLims = 'auto';
        else
            cLims = [-0.245; -0.09];
        end

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end

    disp(' ');
end

% switch format
%     
%     case 'B'
%         
%         if plotInst
%             disp('Presenting Instantaneous Pressure Data...');
%             
%             figHold = fig;
%             
%             for i = 1:nFrames
%                 
%                 if i ~= 1
%                     clf(fig)
%                     fig = figHold;
%                 end
%                 
%                 scalarData = pData.Cp{i};
%                 figTime = num2str(pData.time(i), ['%.', num2str(timePrecision), 'f']);
%                 figName = [caseID, '_Instantaneous_Cp_T', erase(figTime, '.')];
%                 CoM = pData.CoP{i};
%                 figSubtitle = [num2str(pData.time(i), ['%.', num2str(timePrecision), 'f']), ' \it{s}'];
% 
%                 fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
%                                         mapPerim, fig, figName, cMap, geometry, contourlines, ...
%                                         xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
%                                         xLimsPlot, yLimsPlot, zLimsPlot, normalise);
%             end
%             
%         end
%         
% end


%% Local Functions

function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
