%% ParaView Data Initialisation v1.1
% ----
% Initialisation of Exported ParaView Planar Field Data for Further Processing
% ----
% Usage: [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialisePVdata(field, normalise)
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

function [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialisePVdata(field, normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('~/Data/Numerical/OpenFOAM', 'Select Case');

    disp(' ');

    disp(['Case: ', caseFolder]);
    
    % Confirm Support
    if ~contains(caseFolder, ["Run_Test", "Windsor"])
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
                data.(plane).position(:,1) = content(:,1);
                data.(plane).position(:,2) = content(:,2);
                data.(plane).position(:,3) = content(:,3);
                data.(plane).p = content(:,4);

            case 'U'
                data.(plane).position(:,1) = content(:,1);
                data.(plane).position(:,2) = content(:,2);
                data.(plane).position(:,3) = content(:,3);
                data.(plane).u = content(:,4);
                data.(plane).v = content(:,5);
                data.(plane).w = content(:,6);

        end
        
%         % Adjust for Small Values
%         data.(plane).position(abs(data.(plane).position) < 1e-9) = 0;
        
        % Remove Duplicate Entries
        [data.(plane).position, index] = unique(data.(plane).position, 'rows', 'stable');

        switch field

            case 'p'
                data.(plane).p = data.(plane).p(index);

            case 'U'
                data.(plane).u = data.(plane).u(index);
                data.(plane).v = data.(plane).v(index);
                data.(plane).w = data.(plane).w(index);

        end
        
        % Check for Unwanted Planar Deviations
        uniqueX = unique(data.(plane).position(:,1));
        uniqueY = unique(data.(plane).position(:,2));
        uniqueZ = unique(data.(plane).position(:,3));
        
        if height(uniqueX) > 1 && (height(uniqueX) / height(data.(plane).position) < 0.01)
            [count, value] = groupcounts(data.(plane).position(:,1));
            [~, index] = max(count);
            data.(plane).position(:,1) = value(index);
            uniqueX = unique(data.(plane).position(:,1));
        elseif height(uniqueY) > 1 && (height(uniqueY) / height(data.(plane).position) < 0.01)
            [count, value] = groupcounts(data.(plane).position(:,2));
            [~, index] = max(count);
            data.(plane).position(:,2) = value(index);
            uniqueY = unique(data.(plane).position(:,2));
        elseif height(uniqueZ) > 1 && (height(uniqueZ) / height(data.(plane).position) < 0.01)
            [count, value] = groupcounts(data.(plane).position(:,3));
            [~, index] = max(count);
            data.(plane).position(:,3) = value(index);
            uniqueZ = unique(data.(plane).position(:,3));
        end
        
        % Identify Plane Orientation
        if height(uniqueX) == 1
            data.(plane).planeOrientation = 'X';
            data.(plane).planePosition = uniqueX;
        elseif height(uniqueY) == 1
            data.(plane).planeOrientation = 'Y';
            data.(plane).planePosition = uniqueY;
        elseif height(uniqueZ) == 1
            data.(plane).planeOrientation = 'Z';
            data.(plane).planePosition = uniqueZ;
        else
            error('Invalid Dataset (Unsupported Plane Orientation)');
        end
        
    end

    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, precision] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if contains(caseFolder, ["Run_Test", "Windsor"])

            for i = 1:height(dataFiles)
                plane = dataFiles(i).name(1:(end - 4));

                data.(plane).position(:,1) = round(data.(plane).position(:,1) / 1.044, precision);
                data.(plane).position(:,2) = round(data.(plane).position(:,2) / 1.044, precision);
                data.(plane).position(:,3) = round(data.(plane).position(:,3) / 1.044, precision);

                if strcmp(data.(plane).planeOrientation, 'X')
                    data.(plane).planePosition = round(data.(plane).planePosition / 1.044, precision);
                elseif strcmp(data.(plane).planeOrientation, 'Y')
                    data.(plane).planePosition = round(data.(plane).planePosition / 1.044, precision);
                elseif strcmp(data.(plane).planeOrientation, 'Z')
                    data.(plane).planePosition = round(data.(plane).planePosition / 1.044, precision);
                end

            end

        end
        
    end

end