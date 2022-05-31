%% Lagrangian Data Initialisation v1.0
% ----
% Initialisation of OpenFOAM v7 Lagrangian Data for Further Processing
% ----
% Usage: [LagProps, LagDataPlane, LagDataSurface, ...
%           LagDataVolume, sampleInterval] = initialiseLagData(caseFolder, caseName, cloudName, ...
%                                                              plane, surface, volume, ...
%                                                              timeDirs, deltaT, timePrecision, nProc);
%        'caseFolder'    -> Case Path, Stored as s String
%        'caseName'      -> Case Name, Stored as a String
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


%% Main Function

function [LagProps, LagDataPlane, LagDataSurface, ...
          LagDataVolume, sampleInterval] = initialiseLagData(caseFolder, caseName, cloudName, ...
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

    else
        error('Invalid Case Directory (No Volume Data Available)');
    end
    
    if plane
        
        if ~exist([caseFolder, '/LagrangianExtractionPlane'], 'dir')
            error('Invalid Case Directory (No Plane Data Available)');
        end

        dataFilesPlane = dir([caseFolder, '/LagrangianExtractionPlane/LagrangianExtractionPlaneData_*']);

        if isempty(dataFilesPlane)
            error('Invalid Case Directory (No Plane Data Available)');
        end
        
    end

    if surface
        
        if ~exist([caseFolder, '/LagrangianSurfaceContamination'], 'dir')
            error('Invalid Case Directory (No Surface Data Available)');
        end

        dataFilesSurface = dir([caseFolder, '/LagrangianSurfaceContamination/LagrangianSurfaceContaminationData']);

        if isempty(dataFilesSurface)
            error('Invalid Case Directory (No Surface Data Available)');
        end
        
    end
    
    % Initialise Lagrangian Properties
    LagProps = {'d'; 'nParticle'; 'origId'; 'origProcId'; 'positionCartesian'; 'U'};
    LagDataPlane = [];
    LagDataSurface = [];
    LagDataVolume = [];
    
    disp(' ');
    disp(' ');
    
    % Load Previously Collated Data (If Desired/Possible)
    disp('Lagrangian Data Load');
    disp('---------------------');
    
    if plane
        
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Load Saved Plane Data? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/LagData/plane/', caseName, '/*.mat'], ...
                                                 'Select Plane Data');

                if contains(filePath, ['LagData/plane/', caseName])
                    disp(['    Loading ''', fileName, '''...']);
                    LagDataPlane = load([filePath, fileName], 'LagData').LagData;
                    sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
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
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/LagData/surface/', caseName, '/*.mat'], ...
                                                 'Select Surface Data');

                if contains(filePath, ['LagData/surface/', caseName])
                    disp(['    Loading ''', fileName, '''...']);
                    LagDataSurface = load([filePath, fileName], 'LagData').LagData;
                    sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
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
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
                [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/LagData/volume/', caseName, '/*.mat'], ...
                                                 'Select Volume Data');

                if contains(filePath, ['LagData/volume/', caseName])
                    disp(['    Loading ''', fileName, '''...']);
                    LagDataVolume = load([filePath, fileName], 'LagData').LagData;
                    sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
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
    
    if valid
        return
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
                continue
            end
            
            endTime = inputTime('End');
            
            if endTime == -1
                continue
            elseif endTime < startTime
                disp('        WARNING: Invalid Time Format (''endTime'' Precedes ''startTime'')');
                continue
            elseif endTime < str2double(timeDirs(1).name) || startTime > str2double(timeDirs(end).name)
                disp('        WARNING: No Lagrangian Data in Selected Time Range');
                continue
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
                continue
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
    
    % Collate Lagrangian Data
    if plane
        disp(' ');
        LagDataPlane = readLagDataPlane(caseFolder, caseName, LagProps, ...
                                        sampleInterval, timeDirs, deltaT, timePrecision);
    end
    
    if surface
        disp(' ');
        LagDataSurface = readLagDataSurface(caseFolder, caseName, LagProps, ...
                                            sampleInterval, timeDirs, deltaT, timePrecision);
    end
    
    if volume
        disp(' ');
        LagDataVolume = readLagDataVolume(caseFolder, caseName, cloudName, LagProps, ...
                                          sampleInterval, timeDirs, deltaT, timePrecision, nProc);
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
    
    if isnan(newFreq) || length(newFreq) > 1 || newFreq <= 0
        disp('            WARNING: Invalid Entry');
        sampleInterval = -1;
    elseif mod(origFreq, newFreq) ~= 0
        disp(['            WARNING: New Frequency Must Be a Factor of ', num2str(origFreq),' Hz']);
        sampleInterval = -1;
    else
        sampleInterval = origFreq / newFreq;
    end
    
end