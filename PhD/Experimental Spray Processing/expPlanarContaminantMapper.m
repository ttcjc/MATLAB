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

nProc = maxNumCompThreads - 2; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Experimental Spray Mapper v1.1

figSave = false; % Save .fig File(s);

normalise = true;

normalisationValue = 0.031814910694444; % Value Used to Normalise Spray (SB, 1.0 L, 15 Hz)

disp('=====================================');
disp('Planar Experimental Spray Mapper v1.1');
disp('=====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Changes to 'plotPlanarScalarField'


%% Initialise Case

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Initialise Experimental Spray Data
[caseFolder, caseID, expSprayData, samplingFrequency] = initialiseExpSprayData(saveLocation, nProc);

nameDelimiter = strfind(caseID, '_');
planeUnit = strfind(caseID, 'L_');
planeDelimeter = max(nameDelimiter(nameDelimiter < planeUnit));

if contains(caseFolder, 'Far_Field_Soiling_07_22')
    
    if normalise
        planeID = caseID((planeDelimeter + 1):(planeUnit - 1));
    else
        planeID = num2str(str2double(caseID((planeDelimeter + 1):(planeUnit - 1))) * 1.044);
    end
    
end

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
timePrecision = 3;

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

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(expSprayData.time);

mapData.positionGrid = expSprayData.positionGrid;
mapData.time = expSprayData.time;
mapData.inst.seedingDensity = expSprayData.seedingDensity;
clear expSprayData;

% Normalise Position Data
if normalise
    
    if contains(caseFolder, 'Far_Field_Soiling_07_22')
        mapData.positionGrid = round((mapData.positionGrid / 1.044), spacePrecision);
    end
    
end

disp(' ');

% Calculate Time-Averaged Contaminant Map
disp('    Generating Time-Averaged Contaminant Map...');

mapData.mean.seedingDensity = zeros([height(mapData.positionGrid),1]);

for i = 1:nTimes
    mapData.mean.seedingDensity = mapData.mean.seedingDensity + mapData.inst.seedingDensity{i};
end
clear i;

mapData.mean.seedingDensity = mapData.mean.seedingDensity / nTimes;

% Normalise Contaminant Maps
mapData.inst.seedingDensityNorm = mapData.inst.seedingDensity;

for i = 1:nTimes
    mapData.inst.seedingDensityNorm{i} = mapData.inst.seedingDensityNorm{i} / normalisationValue;
end

mapData.mean.seedingDensityNorm = mapData.mean.seedingDensity / normalisationValue;

% Calculate Standard Deviation
mapData.std.seedingDensity = zeros([height(mapData.positionGrid),1]);
mapData.std.seedingDensityNorm = mapData.std.seedingDensity;

for i = 1:nTimes
    mapData.std.seedingDensity = mapData.std.seedingDensity + (mapData.inst.seedingDensity{i} - mapData.mean.seedingDensity).^2;
    mapData.std.seedingDensityNorm = mapData.std.seedingDensityNorm + (mapData.inst.seedingDensityNorm{i} - mapData.mean.seedingDensityNorm).^2;
end

mapData.std.seedingDensity = sqrt((1 / nTimes) * mapData.std.seedingDensity);
mapData.std.seedingDensityNorm = sqrt((1 / nTimes) * mapData.std.seedingDensityNorm);

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

% valid = false;
% while ~valid
%     disp(' ');
%     selection = input('Plot Standard Deviation of Contamination? [y/n]: ', 's');
% 
%     if selection == 'n' | selection == 'N' %#ok<OR2>
%         plotStd = false;
%         valid = true;
%     elseif selection == 'y' | selection == 'Y' %#ok<OR2>
%         plotStd = true;
%         valid = true;
%     else
%         disp('    WARNING: Invalid Entry');
%     end
% 
% end
% clear valid;

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Instantaneous Contamination? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(nTimes);
        
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

if plotMean || plotInst
    
    % Define Plot Limits
    xLimsData = mapData.positionGrid(1,1);
    yLimsData = [min(mapData.positionGrid(:,2)); max(mapData.positionGrid(:,2))];
    zLimsData = [min(mapData.positionGrid(:,3)); max(mapData.positionGrid(:,3))];
    
    if contains(caseFolder, 'Far_Field_Soiling_07_22')
        xLimsPlot = [0.31875; 2.73575];
        yLimsPlot = [-0.522; 0.522];
        zLimsPlot = [0; 0.6264];
%         xLimsPlot = [0.31875; 4.65925];
%         yLimsPlot = [-0.399; 0.218];
%         zLimsPlot = [0.0105; 0.4985];
    end

    if normalise
        
        if contains(caseFolder, 'Far_Field_Soiling_07_22')
            xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
            yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
            zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
        end
        
    end
    
    orientation = 'YZ';
    positionData = mapData.positionGrid;
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
end

% Remove Geometry From Empty Tunnel
if contains(caseID, 'ET')
    geometry = [];
end

if plotMean
    disp('    Presenting Time-Averaged Seeding Density...');

    scalarData = mapData.mean.seedingDensityNorm;
    figName = 'Time_Averaged_Experimental_Seeding_Density';
    contourlines = [0.02; 0.02];
    figSubtitle = ' ';
    cLims = [0; 1.2];

    [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                           mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                           xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
    
    disp(' ');
end

% if plotStd
%     disp('    Presenting Standard Deviation of Seeding Density...');
% 
%     scalarData = mapData.std.seedingDensityNorm;
%     figName = 'Standard_Deviation_Experimental_Seeding_Density';
%     contourlines = [0.012; 0.012];
%     figSubtitle = ' ';
%     cLims = [0; 0.6];
% 
%     [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
%                                            mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
%                                            xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
%                                            xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
%     
%     disp(' ');
% end

if plotInst
    disp('    Presenting Instantaneous Seeding Density...');
    
    contourlines = [0.02; 0.02];
    cLims = [0; 2.2];

    figHold = fig;

    for i = 1:nFrames

        if i ~= 1
            clf(fig);
            fig = figHold;
        end

        scalarData = mapData.inst.seedingDensityNorm{i};
        figTime = num2str(mapData.time(i), '%.2f');
        figName = ['Instantaneous_Experimental_Seeding_Density_T', erase(figTime, '.')];
        figSubtitle = [figTime, ' \it{s}'];
        
        [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
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
        
        if ~exist([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseID], 'dir')
            mkdir([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseID]);
        end
        
        disp(['    Saving to: ', saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseID, '/', dataID '.mat']);
        save([saveLocation, '/Experimental/MATLAB/planarContaminantMap/', caseID, '/', dataID '.mat'], ...
             'caseID', 'planeID', 'dataID', 'mapData', 'sampleInterval', 'timePrecision', 'normalise', '-v7.3', '-noCompression');
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