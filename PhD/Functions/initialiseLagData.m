%% Lagrangian Data Initialisation v1.1
% ----
% Initialisation of OpenFOAM v7 Lagrangian Data for Further Processing
% ----
% Usage: [dataID, LagProps, LagDataPlane, LagDataSurface, ...
%         LagDataVolume, sampleInterval, format] = initialiseLagData(saveLocation, caseFolder, caseID, cloudName, ...
%                                                                    plane, surface, volume, timeDirs, deltaT, ...
%                                                                    timePrecision, nProc);
%        'saveLocation'  -> Start of File Path, Stored as a String
%        'caseFolder'    -> Case Path, Stored as s String
%        'caseID'        -> Case Name, Stored as a String
%        'cloudName'     -> OpenFOAM Cloud Name, Stored as a String
%        'plane'         -> Collect Planar Data [True/False]
%        'surface'       -> Collect Surface Data [True/False]
%        'volume'        -> Collect Volume Data [True/False]
%        'timeDirs'      -> Time Directories, Obtained With 'timeDirectories.m'
%        'deltaT'        -> Time Delta Between Directiories, Obtained With 'timeDirectories.m'
%        'timePrecision' -> Required Rounding Precision for 'deltaT', Obtained With 'timeDirectories.m'
%        'nProc'         -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for 'Snapshot' Data Collation


%% Main Function

function [dataID, LagProps, LagDataPlane, LagDataSurface, ...
          LagDataVolume, sampleInterval, format] = initialiseLagData(saveLocation, caseFolder, caseID, cloudName, ...
                                                                     plane, surface, volume, ...
                                                                     timeDirs, deltaT, timePrecision, nProc)

    % Confirm Lagrangian Data Availability
    i = 1;
    while i <= height(timeDirs)

        if ~exist([caseFolder, '/', timeDirs(i).name, '/lagrangian/', cloudName, '/active'], 'file')
            timeDirs(i) = [];
        else
            i = i + 1;
        end

    end
    clear i;

    if ~isempty(timeDirs)
        disp('Lagrangian Data Identified in the Following Time Directories:');

        for i = 1:height(timeDirs)
            disp(['    ', timeDirs(i).name]);
        end
        clear i;

    else
        error('Invalid Case Directory (No Volume Data Available)');
    end
    
    if plane
        
        if ~exist([caseFolder, '/LagrangianExtractionPlane'], 'dir')
            error('Invalid Case Directory (No Plane Data Available)');
        end

        % Check for Distributed Files
        if ~isempty(dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']))
            distributedFiles = false;
        elseif ~isempty(dir([caseFolder, '/LagrangianExtractionPlane/*/LagrangianExtractionPlaneData_*']))
            distributedFiles = true;
        else
            error('Invalid Case Directory (No Plane Data Available)');
        end
        
    end

    if surface
        
        if ~exist([caseFolder, '/LagrangianSurfaceContamination'], 'dir')
            error('Invalid Case Directory (No Surface Data Available)');
        end

        % Check for Distributed Files
        if ~isempty(dir([caseFolder, '/LagrangianSurfaceContamination/LagrangianSurfaceContaminationData']))
            distributedFiles = false;
        elseif ~isempty(dir([caseFolder, '/LagrangianSurfaceContamination/*/LagrangianSurfaceContaminationData']))
            distributedFiles = true;
        else
            error('Invalid Case Directory (No Surface Data Available)');
        end
        
    end
    
    % Initialise Lagrangian Properties
    LagProps = {'d'; 'nParticle'; 'origId'; 'origProcId'; 'positionCartesian'; 'U'; 'Uslip'};
    LagDataPlane = [];
    LagDataSurface = [];
    LagDataVolume = [];
    
    disp(' ');
    disp(' ');
    
    % Load Previously Collated Data (If Desired/Possible)
    disp('Lagrangian Data Load');
    disp('---------------------');
    
    loadData = true;
    while loadData
        sampleInterval = [];

        if plane
            
            valid = false;
            while ~valid
                disp(' ');
                selection = input('Load Saved Plane Data? [y/n]: ', 's');
    
                if selection == 'n' | selection == 'N' %#ok<OR2>
                    break;
                elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/plane/*.mat'], ...
                                                     'Select Plane Data');
    
                    if contains(filePath, ['/LagData/', caseID, '/plane'])
                        disp(['    Loading ''', fileName, '''...']);
                        
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        LagDataPlane = load([filePath, fileName], 'LagData').LagData;
                        sampleInterval = [sampleInterval; load([filePath, fileName], 'sampleInterval').sampleInterval]; %#ok<AGROW>
                        formatPlane = load([filePath, fileName], 'format').format;
                        
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
            clear valid;
            
        end
        
        if surface
            
            valid = false;
            while ~valid
                disp(' ');
                selection = input('Load Saved Surface Data? [y/n]: ', 's');
    
                if selection == 'n' | selection == 'N' %#ok<OR2>
                    break;
                elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/surface/*.mat'], ...
                                                     'Select Surface Data');
    
                    if contains(filePath, ['/LagData/', caseID, '/surface'])
                        disp(['    Loading ''', fileName, '''...']);
                        
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        LagDataSurface = load([filePath, fileName], 'LagData').LagData;
                        sampleInterval = [sampleInterval; load([filePath, fileName], 'sampleInterval').sampleInterval]; %#ok<AGROW>
                        formatSurface = load([filePath, fileName], 'format').format;
                        
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
            clear valid;
            
        end
        
        if volume
            
            valid = false;
            while ~valid
                disp(' ');
                selection = input('Load Saved Volume Data? [y/n]: ', 's');
    
                if selection == 'n' | selection == 'N' %#ok<OR2>
                    break;
                elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/LagData/', caseID, '/volume/*.mat'], ...
                                                     'Select Volume Data');
    
                    if contains(filePath, ['/LagData/', caseID, '/volume'])
                        disp(['    Loading ''', fileName, '''...']);
                        
                        dataID = load([filePath, fileName], 'dataID').dataID;
                        LagDataVolume = load([filePath, fileName], 'LagData').LagData;
                        sampleInterval = [sampleInterval; load([filePath, fileName], 'sampleInterval').sampleInterval]; %#ok<AGROW>
                        
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
            clear valid;
            
        end
        
        loadData = false;
    end
    clear loadData;
    
    % Confirm Required Data Exists
    valid = true;

    if plane && isempty(LagDataPlane)
        valid = false;
    end
    
    if surface && isempty(LagDataSurface)
        valid = false;
    end
    
    if volume && isempty(LagDataVolume)
        valid = false;
    end

    % Confirm Matching Sampling Intervals
    if valid
        
        if height(unique(sampleInterval)) ~= 1
            valid = false;
        else
            sampleInterval = sampleInterval(1);
        end
        
    end

    % Confirm Matching Formats
    if valid && plane && surface
        valid = strcmp(formatPlane, formatSurface);
    end

    % Return if Data Required Data Is Present and Valid
    if valid
        return;
    elseif height(unique(sampleInterval)) ~= 1
        disp(' ');
        
        disp('WARNING: Loaded Sample Intervals Do Not Match');
        disp('         Collating New Data...');
    elseif plane && surface 
        disp(' ');
        
        disp('WARNING: Loaded Plane and Surface Datasets Have Different Formats');
        disp('         Collating New Data...');
    else
        disp(' ');
        
        disp('WARNING: Required Data Unavailable');
        disp('         Collating New Data...');
    end
    
    disp(' ');
    disp(' ');
    
    % Acquire Lagrangian Data
    disp('Lagrangian Data Acquisition');
    disp('----------------------------');
    
    disp(' ');

    disp('Storing the Following Lagrangian Properties:')
    disp('    d');
    disp('    nParticle');
    disp('    origId');
    disp('    origProcId');
    disp('    positionCartesian');
    disp('    U');
    disp('    Uslip');
    
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
            elseif endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Lagrangian Data in Selected Time Range');
                continue;
            end

            i = 1;
            while i <= height(timeDirs)
                
                if str2double(timeDirs(i).name) < startTime || str2double(timeDirs(i).name) > endTime
                    timeDirs(i) = [];
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
        disp(['Default Sampling Frequency: ', num2str(round((1 / deltaT), timePrecision)), ' Hz']);
        selection = input('    Reduce Recording Frequency? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            sampleInterval = 1;
            
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            sampleInterval = inputFreq(round((1 / deltaT), timePrecision));
            
            if sampleInterval == -1
                continue;
            end
            
            if sampleInterval >= (floor(height(timeDirs) / 2))
                disp('            WARNING: Sampling Interval Must Fall Within Data Range');
            else
                valid = true;
            end
            
        else
            disp('        WARNING: Invalid Entry');
        end

    end
    clear valid;

    % Specify Data Collation Format
    if plane || surface
        disp(' ');
        
        disp('Possible Data Collation Formats:');
        disp('    A: Cumulative');
        disp('    B: Snapshot');
        
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Select Data Collation Format [A/B]: ', 's');
            
            if selection == 'a' | selection == 'A' %#ok<OR2>
                format = 'cumulative';
                valid = true;
            elseif selection == 'b' | selection == 'B' %#ok<OR2>
                format = 'snapshot';
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end
        
        end
        clear valid;

    end
    
    % Define Data ID
    startInst = erase(num2str(str2double(timeDirs(1).name), ['%.', num2str(timePrecision), 'f']), '.');
    endInst = erase(num2str(str2double(timeDirs(end).name), ['%.', num2str(timePrecision), 'f']), '.');
    freq = num2str(round((1 / (deltaT * sampleInterval)), timePrecision));

    dataID = ['T', startInst, '_T', endInst, '_F', freq];
    
    % Collate Lagrangian Data
    if plane
        disp(' ');
        
        LagDataPlane = readLagDataPlane(saveLocation, caseFolder, caseID, distributedFiles, dataID, ...
                                        LagProps, timeDirs, sampleInterval, format);
    end
    
    if surface
        disp(' ');
        
        LagDataSurface = readLagDataSurface(saveLocation, caseFolder, caseID, distributedFiles, dataID, ...
                                            LagProps, timeDirs, sampleInterval, format);
    end
    
    if volume
        disp(' ');
        
        LagDataVolume = readLagDataVolume(saveLocation, caseFolder, caseID, dataID, cloudName, LagProps, ...
                                          timeDirs, sampleInterval, nProc);
    end
    
    % Update Data ID
    if plane || surface
        dataID = [dataID, '_', format];
    end

end


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