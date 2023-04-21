%% Planar Experimental Spray Mapper v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

normalise = false; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

normalisationValue = 0.031814910694444; % SB_1.0L_120s_15Hz Time-Averaged Max

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('=====================================');
disp('Planar Experimental Spray Mapper v1.0');
disp('=====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Case

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');

% Load Experimental Spray Data
[caseFolder, caseName, expSprayData, samplingFrequency] = initialiseExpSprayData(saveLocation, nProc);

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

% Select Times of Interest
valid = false;
while ~valid
    disp(' ');
    selection = input('Restrict Data Range? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        startTime = inputTime('Start');

        if startTime == -1
            continue;
        end

        endTime = inputTime('End');

        if endTime == -1
            continue;
        elseif endTime < startTime
            disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
            continue;
        elseif endTime < expSprayData.time(1) || startTime > expSprayData.time(end)
            disp('        WARNING: No Lagrangian Data in Selected Time Range');
            continue;
        end

        i = 1;
        while i <= height(expSprayData.time)

            if expSprayData.time(i) < startTime || expSprayData.time(i) > endTime
                expSprayData.time(i) = [];
                expSprayData.seedingDensity(i) = [];
            else
                i = i + 1;
            end

        end
        clear i;

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

% Specify Sampling Frequency
valid = false;
while ~valid
    disp(' ');
    disp(['Default Sampling Frequency: ', num2str(samplingFrequency), ' Hz']);
    selection = input('    Reduce Recording Frequency? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        sampleInterval = 1;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        sampleInterval = inputFreq(samplingFrequency);

        if sampleInterval == -1
            continue;
        end

        if sampleInterval >= (floor(height(expSprayData.time) / 2))
            disp('            WARNING: Sampling Interval Must Fall Within Data Range');
        else
            valid = true;
        end

    else
        disp('        WARNING: Invalid Entry');
    end

end
clear valid;

if sampleInterval ~= 1
    expSprayData.time = expSprayData.time(1:sampleInterval:end);
    expSprayData.seedingDensity = expSprayData.seedingDensity(1:sampleInterval:end);
end

% Define Data ID
timePrecision = 2;

startInst = erase(num2str(expSprayData.time(1), ['%.', num2str(timePrecision), 'f']), '.');
endInst = erase(num2str(expSprayData.time(end), ['%.', num2str(timePrecision), 'f']), '.');

samplingFrequency = num2str(samplingFrequency / sampleInterval);

if normalise
    dataID = ['T', startInst, '_T', endInst, '_F', samplingFrequency, '_Norm'];
else
    dataID = ['T', startInst, '_T', endInst, '_F', samplingFrequency];
end

disp(' ');
disp(' ');


%% Generate Contaminant Maps

disp('Contaminant Mapping');
disp('--------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

mapData.positionGrid = expSprayData.positionGrid;
mapData.inst.time = expSprayData.time;
mapData.inst.seedingDensity = expSprayData.seedingDensity;
clear expSprayData;

% ---- Temporary ---- %
[mapData.positionGrid, index] = sortrows(mapData.positionGrid,3);

for i = 1:height(mapData.inst.time)
    mapData.inst.seedingDensity{i} = mapData.inst.seedingDensity{i}(index);
end
% ---- Temporary ---- %

% Normalise Position Data
if normalise
    mapData.positionGrid = round((mapData.positionGrid / 1.044), spacePrecision);
end

disp(' ');

% % Calculate Instantaneous Centre of Mass
% disp('    Calculating Instantaneous Centre of Mass...');
% 
% % Initialise Progress Bar
% wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
% wB.Children.Title.Interpreter = 'none';
% 
% mapData.inst.CoM = cell(height(mapData.inst.time),1);
% 
% for i = 1:height(mapData.inst.time)
%     mapData.inst.CoM{i} = zeros(1,3);
%     
%     mapData.inst.CoM{i}(1) = mapData.positionGrid(1,1);
%     mapData.inst.CoM{i}(2) = sum(mapData.inst.seedingDensity{i} .* mapData.positionGrid(:,2)) / sum(mapData.inst.seedingDensity{i});
%     mapData.inst.CoM{i}(3) = sum(mapData.inst.seedingDensity{i} .* mapData.positionGrid(:,3)) / sum(mapData.inst.seedingDensity{i});
%     
%     waitbar((i / height(mapData.inst.time)), wB);
% end
% 
% delete(wB);
% 
% disp(' ');

% Calculate Time-Averaged Contaminant Map
disp('    Generating Time-Averaged Contaminant Map...');

mapData.mean.seedingDensity = zeros(height(mapData.positionGrid),1);

for i = 1:height(mapData.inst.time)
    mapData.mean.seedingDensity = mapData.mean.seedingDensity + mapData.inst.seedingDensity{i};
end

mapData.mean.seedingDensity = mapData.mean.seedingDensity / height(mapData.inst.time);

% Normalise Contaminant Maps
if normalise

    for i = 1:height(mapData.inst.time)
        mapData.inst.seedingDensity{i} = mapData.inst.seedingDensity{i} / normalisationValue;
    end
    
    mapData.mean.seedingDensity = mapData.mean.seedingDensity / normalisationValue;
end

% % Calculate Time-Averaged Centre of Mass
% mapData.mean.CoM = zeros(1,3);
%     
% mapData.mean.CoM(1) = mapData.positionGrid(1,1);
% mapData.mean.CoM(2) = sum(mapData.mean.seedingDensity .* mapData.positionGrid(:,2)) / sum(mapData.mean.seedingDensity);
% mapData.mean.CoM(3) = sum(mapData.mean.seedingDensity .* mapData.positionGrid(:,3)) / sum(mapData.mean.seedingDensity);

% Calculate Standard Deviation
mapData.std.seedingDensity = zeros(height(mapData.positionGrid),1);

for i = 1:height(mapData.inst.time)
    mapData.std.seedingDensity = mapData.std.seedingDensity + (mapData.inst.seedingDensity{i} - mapData.mean.seedingDensity).^2;
end

mapData.std.seedingDensity = sqrt((1 / height(mapData.inst.time)) * mapData.std.seedingDensity);

evalc('delete(gcp(''nocreate''));');
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
    selection = input('Plot Time-Averaged Contamination? [y/n]: ', 's');

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

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Standard Deviation of Contamination? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotStd = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotStd = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Instantaneous Contamination? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(height(mapData.inst.time));
        
        if nFrames == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotStd || plotInst    
    % Define Plot Limits
    xLimsData = mapData.positionGrid(1,1);
    yLimsData = [min(mapData.positionGrid(:,2)); max(mapData.positionGrid(:,2))];
    zLimsData = [min(mapData.positionGrid(:,3)); max(mapData.positionGrid(:,3))];
    
%     % Plot Numerical Data Range
%     xLimsPlot = [0.31875; 4.65925];
%     yLimsPlot = [-0.5945; 0.5945];
%     zLimsPlot = [0; 0.739];
    
    % Plot Experimental Data Range
    xLimsPlot = [0.31875; 4.65925];
    yLimsPlot = [-0.399; 0.218];
    zLimsPlot = [0.0105; 0.4985];

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end
    
    orientation = 'YZ';
    positionData = mapData.positionGrid;
    mapPerim = [];
    cMap = flipud(viridis(32));
    CoM = [];
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    nPlanes = 1;
    planeNo = 1;
end

% Remove Geometry From Empty Tunnel
if contains(caseName, 'ET')
    geometry = [];
end

if plotMean
    disp('    Presenting Time-Averaged Contamination...');

    scalarData = mapData.mean.seedingDensity;
    figName = [caseName, '_Experimental_Contamination_Time_Average'];
%     contourlines = (0.2:0.2:0.8);
    contourlines = [];
    figSubtitle = ' ';
    cLims = [0; 0.032];

    fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                            mapPerim, fig, figName, cMap, geometry, contourlines, ...
                            xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                            xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
    
    disp(' ');
end

if plotStd
    disp('    Presenting Standard Deviation of Contamination...');

    scalarData = mapData.std.seedingDensity;
    figName = [caseName, '_Experimental_Contamination_Standard_Deviation'];
    contourlines = [];
    figSubtitle = ' ';
    cLims = [0; 0.6];

    fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                            mapPerim, fig, figName, cMap, geometry, contourlines, ...
                            xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                            xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
    
    disp(' ');
end

if plotInst
    disp('    Presenting Instantaneous Contamination...');
    
    contourlines = [];
    cLims = [0; 0.1];
%     cLims = [0; max(cellfun(@max, mapData.inst.seedingDensity))];

    figHold = fig;

    for i = 1:nFrames

        if i ~= 1
            clf(fig);
            fig = figHold;
        end

        scalarData = mapData.inst.seedingDensity{i};
        figTime = num2str(mapData.inst.time(i), '%.2f');
        figName = [caseName, '_Experimental_Contamination_Instantaneous_T', erase(figTime, '.')];
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
    end
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Map Presentation');
    
    disp(' ');
end

disp(' ');


%% Save Map Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        if ~exist([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseName], 'dir')
            mkdir([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseName]);
        end
        
        disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseName, '/', dataID '.mat']);
        save([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseName, '/', dataID '.mat'], ...
             'dataID', 'mapData', 'sampleInterval', 'normalise', 'timePrecision', '-v7.3', '-noCompression');
        disp('        Success');
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function time = inputTime(type)

    time = str2double(input(['    Input ', type, ' Time [s]: '], 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        time = -1;
    end

end


function sampleInterval = inputFreq(origFreq)
    
    newFreq = str2double(input('        Input Frequency [Hz]: ', 's'));
    
    if isnan(newFreq) || newFreq <= 0 || newFreq > origFreq
        disp('            WARNING: Invalid Entry');
        sampleInterval = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        sampleInterval = -1;
    else
        sampleInterval = origFreq / newFreq;
    end
    
end


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end
