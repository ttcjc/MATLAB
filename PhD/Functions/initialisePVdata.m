%% ParaView Data Initialisation v1.1
% ----
% Initialisation of Exported ParaView Planar Field Data for Further Processing
% ----
% Usage: [caseFolder, PVdata, geometry, ...
%         xDims, yDims, zDims, spacePrecision] = initialisePVdata(field, normalise)
%        'field'     -> Desired Field Stored as String
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [caseName, PVdata, geometry, ...
          xDims, yDims, zDims, spacePrecision] = initialisePVdata(field, normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('/mnt/Processing/Data/Numerical/OpenFOAM', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseName = caseFolder(namePos:end);

    disp(' ');

    disp(['Case: ', caseName]);
    
    % Confirm Support
    if ~contains(caseName, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    % Confirm Data Availability
    if ~strcmp(field, 'p') && ~strcmp(field, 'U')
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end

    switch field

        case 'p'
            fieldLabel = 'Pressure';
            dataFiles = dir([caseFolder, '/p_*.csv']);

        case 'U'
            fieldLabel = 'Velocity';
            dataFiles = dir([caseFolder, '/U_*.csv']);

    end
    
    if isempty(dataFiles)
        error(['Invalid Case Directory (No ', fieldLabel, ' Data Found)']);
    end
    
    disp(' ');
    disp(' ');

    % Acquire ParaView Data
    disp('ParaView Data Acquisition');
    disp('--------------------------');
    
    disp(' ');

    disp('Formatting:');

    for i = 1:height(dataFiles)
        plane = dataFiles(i).name(1:(end - 4)); % Ignore .csv
        disp(['    ', plane]);
    
        % Load Data
        content = readmatrix([caseFolder, '/', dataFiles(i).name]);
        
        % Format Data
        switch field

            case 'p'
                PVdata.(plane).position = content(:,[1,2,3]);
                PVdata.(plane).pMean = content(:,4);

            case 'U'
                PVdata.(plane).position = content(:,[1,2,3]);
                PVdata.(plane).uMean = content(:,4);
                PVdata.(plane).vMean = content(:,5);
                PVdata.(plane).wMean = content(:,6);

        end
        
        % Check for Unwanted Planar Deviations
        uniqueX = unique(PVdata.(plane).position(:,1));
        uniqueY = unique(PVdata.(plane).position(:,2));
        uniqueZ = unique(PVdata.(plane).position(:,3));
        
        if height(uniqueX) > 1 && (height(uniqueX) / height(PVdata.(plane).position) < 0.01)
            [count, value] = groupcounts(PVdata.(plane).position(:,1));
            [~, index] = max(count);
            PVdata.(plane).position(:,1) = value(index);
            uniqueX = unique(PVdata.(plane).position(:,1));
        elseif height(uniqueY) > 1 && (height(uniqueY) / height(PVdata.(plane).position) < 0.01)
            [count, value] = groupcounts(PVdata.(plane).position(:,2));
            [~, index] = max(count);
            PVdata.(plane).position(:,2) = value(index);
            uniqueY = unique(PVdata.(plane).position(:,2));
        elseif height(uniqueZ) > 1 && (height(uniqueZ) / height(PVdata.(plane).position) < 0.01)
            [count, value] = groupcounts(PVdata.(plane).position(:,3));
            [~, index] = max(count);
            PVdata.(plane).position(:,3) = value(index);
            uniqueZ = unique(PVdata.(plane).position(:,3));
        end
        
        % Identify Plane Orientation
        if height(uniqueX) == 1
            PVdata.(plane).planeOrientation = 'YZ';
            PVdata.(plane).planePosition = uniqueX;
        elseif height(uniqueY) == 1
            PVdata.(plane).planeOrientation = 'XZ';
            PVdata.(plane).planePosition = uniqueY;
        elseif height(uniqueZ) == 1
            PVdata.(plane).planeOrientation = 'XY';
            PVdata.(plane).planePosition = uniqueZ;
        else
            error('Invalid Dataset (Unsupported Plane Orientation)');
        end
        
    end

    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if contains(caseName, ["Run_Test", "Windsor"])

            for i = 1:height(dataFiles)
                plane = dataFiles(i).name(1:(end - 4));

                PVdata.(plane).position = round((PVdata.(plane).position / 1.044), spacePrecision);
                PVdata.(plane).planePosition = round((PVdata.(plane).planePosition / 1.044), spacePrecision);
            end

        end
        
    end
    
    % Remove Duplicate Entries
    for i = 1:height(dataFiles)
        plane = dataFiles(i).name(1:(end - 4));
        
        [PVdata.(plane).position, index] = unique(PVdata.(plane).position, 'rows', 'stable');
        
        switch field
            
            case 'p'
                PVdata.(plane).pMean = PVdata.(plane).pMean(index);

            case 'U'
                PVdata.(plane).uMean = PVdata.(plane).uMean(index);
                PVdata.(plane).vMean = PVdata.(plane).vMean(index);
                PVdata.(plane).wMean = PVdata.(plane).wMean(index);
        
        end
        
        PVdata.(plane) = orderfields(PVdata.(plane));
    end
    
    PVdata = orderfields(PVdata);

end