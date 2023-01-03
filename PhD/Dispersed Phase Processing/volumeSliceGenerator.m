%% Lagrangian Volume Slice Generator v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

% saveLocation = '/mnt/Processing/Data';
saveLocation = '~/Data';

multiSlice = true; % Allow Selection of Multiple Slices

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('===========================');
disp('Volume Slice Generator v1.0');
disp('===========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near-Field');
disp('    B: Far-Field');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Region of Interest [A/B]: ', 's');

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
clear valid;

disp(' ');
disp(' ');


%% Acquire Volume Field

disp('Volume Field Acquisition');
disp('-------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/volumeField/*.mat'], ...
                                      'Select Volumetric Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/nearField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/farField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
end
clear valid;

namePos = strfind(filePath, '/');
caseName = filePath((namePos(end - 2) + 1):(namePos(end - 1) - 1));

timePrecision = strfind(fileName, '_T') - 3;

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Identify Volume Slice(s)

% Restore Original Dimensions
if normalise
    volumeData.positionGrid = round((volumeData.positionGrid * 1.044), spacePrecision);
end

% Identify Desired Slice(s)
volumeSlice = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, normalise, multiSlice);

slices = fieldnames(volumeSlice);

% Normalise Dimensions
if normalise
    volumeData.positionGrid = round((volumeData.positionGrid / 1.044), spacePrecision);

    for i = 1:height(slices)
        volumeSlice.(slices{i}).xLims = round((volumeSlice.(slices{i}).xLims / 1.044), spacePrecision);
        volumeSlice.(slices{i}).yLims = round((volumeSlice.(slices{i}).yLims / 1.044), spacePrecision);
        volumeSlice.(slices{i}).zLims = round((volumeSlice.(slices{i}).zLims / 1.044), spacePrecision);
        volumeSlice.(slices{i}).position = round((volumeSlice.(slices{i}).position / 1.044), spacePrecision);
    end

end

disp(' ');
disp(' ');


%% Extract Slice Data

disp('Volume Slice Data Extraction');
disp('-----------------------------');

disp(' ');

for i = 1:height(slices)

    switch volumeSlice.(slices{i}).orientation

        case 'YZ'
            index = find((volumeData.positionGrid(:,1) == volumeSlice.(slices{i}).xLims) & ...
                         (volumeData.positionGrid(:,2) >= volumeSlice.(slices{i}).yLims(1) & volumeData.positionGrid(:,2) <= volumeSlice.(slices{i}).yLims(2)) & ...
                         (volumeData.positionGrid(:,3) >= volumeSlice.(slices{i}).zLims(1) & volumeData.positionGrid(:,3) <= volumeSlice.(slices{i}).zLims(2)));

        case 'XZ'
            index = find((volumeData.positionGrid(:,1) >= volumeSlice.(slices{i}).xLims(1) & volumeData.positionGrid(:,1) <= volumeSlice.(slices{i}).xLims(2)) & ...
                             (volumeData.positionGrid(:,2) == volumeSlice.(slices{i}).yLims) & ...
                             (volumeData.positionGrid(:,3) >= volumeSlice.(slices{i}).zLims(1) & volumeData.positionGrid(:,3) <= volumeSlice.(slices{i}).zLims(2)));

        case 'XY'
            index = find((volumeData.positionGrid(:,1) >= volumeSlice.(slices{i}).xLims(1) & volumeData.positionGrid(:,1) <= volumeSlice.(slices{i}).xLims(2)) & ...
                             (volumeData.positionGrid(:,2) >= volumeSlice.(slices{i}).yLims(1) & volumeData.positionGrid(:,2) <= volumeSlice.(slices{i}).yLims(2)) & ...
                             (volumeData.positionGrid(:,3) == volumeSlice.(slices{i}).zLims));

    end

    % Collate Slice Data
    volumeSlice.(slices{i}).positionGrid = volumeData.positionGrid(index,:);
    volumeSlice.(slices{i}).inst.time = volumeData.inst.time;
    
    for j = 1:height(volumeSlice.(slices{i}).inst.time)
        volumeSlice.(slices{i}).inst.nParticles{j} = volumeData.inst.nParticles{j}(index);
        volumeSlice.(slices{i}).inst.volFraction{j} = volumeData.inst.volFraction{j}(index);
        volumeSlice.(slices{i}).inst.mass{j} = volumeData.inst.mass{j}(index);
%         volumeSlice.(slices{i}).inst.d43{j} = volumeData.inst.d43{j}(index);
%         volumeSlice.(slices{i}).inst.d32{j} = volumeData.inst.d32{j}(index);
%         volumeSlice.(slices{i}).inst.d30{j} = volumeData.inst.d30{j}(index);
%         volumeSlice.(slices{i}).inst.d20{j} = volumeData.inst.d20{j}(index);
%         volumeSlice.(slices{i}).inst.d10{j} = volumeData.inst.d10{j}(index);
    end

        volumeSlice.(slices{i}).mean.nParticles = volumeData.mean.nParticles(index);
        volumeSlice.(slices{i}).mean.volFraction = volumeData.mean.volFraction(index);
        volumeSlice.(slices{i}).mean.mass = volumeData.mean.mass(index);
%         volumeSlice.(slices{i}).mean.d43 = volumeData.mean.d43(index);
%         volumeSlice.(slices{i}).mean.d32 = volumeData.mean.d32(index);
%         volumeSlice.(slices{i}).mean.d30 = volumeData.mean.d30(index);
%         volumeSlice.(slices{i}).mean.d20 = volumeData.mean.d20(index);
%         volumeSlice.(slices{i}).mean.d10 = volumeData.mean.d10(index);
end


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

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        nFrames = inputFrames(height(volumeSlice.(slices{1}).inst.time));
        
        if nFrames == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotInst || plotMean

    if height(slices) > 1
    
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Plot Planes in Single Figure? [y/n]: ', 's');
        
            if selection == 'n' | selection == 'N' %#ok<OR2>
                multiPlane = false;
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                multiPlane = true;
                
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end
        
        end
        clear valid;
    
    end
    
    % Select Variable(s) of Interest
    plotVars = fieldnames(volumeSlice.(slices{1}).mean);

    valid = false;
    while ~valid
        disp(' ');
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', plotVars);

        if ~valid
            disp('WARNING: No Mapping Variables Selected');
        end
    end
    clear valid;

    plotVars = plotVars(index);
    
    disp(['Plotting ', num2str(height(plotVars)), ' Variable(s) of Interest']);
end

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotInst || plotMean
    cMap = flipud(viridis(32));
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    mapPerim = [];
    contourlines = [];
    CoM = [];
    figSubtitle = ' ';

    % Define Plot Limits
    switch format

        case 'A'
            
            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.26625];
                yLimsPlot = [-0.3445; 0.3445];
                zLimsPlot = [0; 0.489];
            end

        case 'B'
            
            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 2.73575]; % 2L + Overhang
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

end

if plotMean

    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);

        planeNo = 1;
        cLims = [0;0];

        if multiPlane
            figName = ['Time_Averaged_', strjoin(slices, '_'), '_', plotVars{i}, '_Map', '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
            nPlanes = height(slices);
        else
            nPlanes = 1;
        end

        for j = 1:height(slices)
            orientation = volumeSlice.(slices{j}).orientation;
            positionData = volumeSlice.(slices{j}).positionGrid;
            xLimsData = volumeSlice.(slices{j}).xLims;
            yLimsData = volumeSlice.(slices{j}).yLims;
            zLimsData = volumeSlice.(slices{j}).zLims;
            scalarData = volumeSlice.(slices{j}).mean.(plotVars{i});

            if ~multiPlane
                figName = ['Time_Averaged_', slices{j}, '_', plotVars{i}, '_Map', '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
                cLims = [0; max(scalarData)];
            else
                cLims = [0; max(max(cLims), max(scalarData))];
            end
            
            [fig, planeNo] = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                               mapPerim, fig, figName, cMap, geometry, contourlines, ...
                                               xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
                                               xLimsPlot, yLimsPlot, zLimsPlot, normalise, nPlanes, planeNo);
        end

    end

end

if ~plotMean && ~plotInst
    disp('    Skipping Volume Field Presentation');

    disp(' ');
end

disp(' ');


















% if plotInst
%     
%     for i = 1:height(plotVars)
%         disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
%         
%         contourlines = [];
%         
%         if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
%             cLims = dLims;
%         elseif strcmp(plotVars{i}, 'massNorm')
%             cLims = [0; 3.6];
%         else
%             cLims = [0; max(cellfun(@max, mapData.inst.(plotVars{i})))];
%         end
%         
%         figHold = fig;
%         
%         for j = 1:nFrames
%             
%             if j ~= 1
%                 clf(fig);
%                 fig = figHold;
%             end
%             
%             scalarData = mapData.inst.(plotVars{i}){j};
%             figTime = num2str(mapData.inst.time(j), ['%.', num2str(timePrecision), 'f']);
%             
%             switch format
%                 
%                 case 'A'
%                     figName = ['Instantaneous_Base_', plotVars{i}, '_Map_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
%                 
%                 case 'B'
%                     figName = ['Instantaneous_', planePos, '_', plotVars{i}, '_Map_T', erase(figTime, '.'), '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
%             end
%             
%         if contains(plotVars{i}, ["mass", "massNorm"])
%             CoM = mapData.inst.CoM{j};
%         else
%             CoM = [];
%         end
%             
%         figSubtitle = [figTime, ' \it{s}'];
%         
%         fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
%                                 mapPerim, fig, figName, cMap, geometry, contourlines, ...
%                                 xDims, yDims, zDims, CoM, figTitle, figSubtitle, cLims, ...
%                                 xLimsPlot, yLimsPlot, zLimsPlot, normalise);
%         end
%         
%     end
%     
%     disp(' ');
% end


%% Save Map Data

