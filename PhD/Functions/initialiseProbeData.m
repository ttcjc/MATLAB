%% Probe Data Initialisation v2.0
% ----
% Initialisation of OpenFOAM v7 Probe Data for Further Processing
% ----
% Usage: [] = initialiseProbeData();

%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v1.2 - Added Support for Volumetric or Planar Data Extraction
% v2.0 - Update To Merge Functionality of 'initialiseVelocityProbeData' and 'initialisePressureProbeData'


%% Supported OpenFOAM Campaigns

% Windsor_Upstream_2023
% Windsor_fullScale


%% Supported Processing Formats

% Volumetric:       'volume'
% Multiple Planes:  'multiPlane'
% Isolated Plane:   'singlePlane'


%% Main Function

%#ok<*OR2>

function [caseFolder, campaignID, caseID, timeDirs, deltaT, ...
          timePrecision, dataID, probeData, sampleInt] = initialiseProbeData(saveLoc, nProc, probeType, ...
                                                                             probeRegion, extractFormat)
    
    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');
    
    campaignID = caseFolder((strfind(caseFolder, 'results/') + 8):(max(strfind(caseFolder, '/')) - 1));
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', campaignID, ', ', caseID]);

    % Confirm Support
    if ~strcmp(campaignID, 'Windsor_Upstream_2023') && ~strcmp(campaignID, 'Windsor_fullScale')
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    disp(' ');
    
    % Confirm Basic Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, probeType);
    
    disp(' ');
    disp(' ');

    % Acquire Probe Data
    disp('Probe Data Acquisition');
    disp('-----------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        
        selection = input('Load Saved Probe Data? [y/n]: ', 's');
    
        if selection == 'n' | selection == 'N'
            disp(' ');
            
            [dataID, probeData, sampleInt] = readProbeData(saveLoc, caseFolder, campaignID, caseID, ...
                                                           timeDirs, deltaT, timePrecision, ...
                                                           probeType, probeRegion, nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y'
            [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/probeData/', campaignID, '/', ...
                                             caseID, '/', probeType, '/', probeRegion, '/*.mat'], ...
                                             'Select Probe Data');
    
            if contains(filePath, caseID) && contains(filePath, probeType) && ...
               contains(filePath, probeRegion)
                disp(['    Loading ''', fileName, '''...']);
                
                dataID = load([filePath, fileName], 'dataID').dataID;
                probeData = load([filePath, fileName], 'probeData').probeData;
                sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
                
                disp('        Success');
                
                valid = true;
            else
                disp('    WARNING: Invalid File Selection');
                clear fileName filePath;
            end
    
        else
            disp('    WARNING: Invalid Entry');
            clear fileName filePath;
        end
    
    end
    
    if strcmp(extractFormat, 'volume')
        return;
    end
    
    disp(' ');
    disp(' ');
    
    % Extract Planar Probe Data
    if strcmp(extractFormat, 'singlePlane')
        volumeSlice = identifyVolumeSlices(probeData.positionGrid, spacePrecision, false, false);
    else
        volumeSlice = identifyVolumeSlices(probeData.positionGrid, spacePrecision, true, false);
    end
    
    switch orientation
        
        case 'YZ'
            index = find((probeData.positionGrid(:,1) == volumeSlice.xLims) & ...
                         ((probeData.positionGrid(:,2) >= volumeSlice.yLims(1)) & ...
                          (probeData.positionGrid(:,2) <= volumeSlice.yLims(2))) & ...
                         ((probeData.positionGrid(:,3) >= volumeSlice.zLims(1)) & ...
                          (probeData.positionGrid(:,3) <= volumeSlice.zLims(2))));
        
        case 'XZ'
        
        case 'XY'
            
    end
    
    % Collate Planar Probe Data
    probeData.positionGrid = probeData.positionGrid(index,:);
    probeData.u = probeData.u(index);
    probeData.v = probeData.v(index);
    probeData.w = probeData.w(index);
    
end