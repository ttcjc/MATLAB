%% ParaView Data Initialisation v2.0
% ----
% Initialisation of Exported ParaView Planar Field Data for Further Processing
% ----
% Usage: [caseID, PVdata, geometry, ...
%         xDims, yDims, zDims, spacePrecision, ...
%         normDims, normLength] = initialisePVdata(field, normDims)
%        'saveLocation' -> Start of File Path, Stored as a String
%        'field'        -> Desired Field Stored as String
%        'normDims'     -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v2.0 - Update To Improve Consistency of Structures Across Repository


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022
% Windsor_Upstream_2023
% Windsor_fullscale


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [campaignID, caseID, PVdata] = initialisePVdata(saveLocation, field)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir([saveLocation, '/Numerical/ParaView'], 'Select Case');
    
    campaignID = caseFolder((strfind(caseFolder, 'ParaView/') + 9):(max(strfind(caseFolder, '/')) - 1));
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', caseID]);
    
    % Confirm Support
    if ~contains(caseID, ["Run_Test", "Windsor"])
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
    
    % Confirm Valid Data Orientation
    for i = 1:height(dataFiles)
        plane = dataFiles(i).name(1:(end - 4)); % Ignore .csv
        
        if contains(plane, 'Base') || contains(plane, '_X_')
            orientation = 'YZ';
        elseif contains(plane, '_Y_')
            orientation = 'XZ';
        elseif contains(plane, '_Z_')
            orientation = 'XY';
        else
            error(['''', plane, ''' Is an Invalid Data File (Unexpected Naming Convention)']);
        end
        
    end
            
    
    disp(' ');
    disp(' ');

    % Acquire ParaView Data
    disp('ParaView Data Acquisition');
    disp('--------------------------');
    
    disp(' ');

    disp('Formatting:');

    for i = 1:height(dataFiles)
        plane = dataFiles(i).name(1:(end - 4));
        
        disp(['    ', plane]);
    
        % Load Data
        content = readmatrix([caseFolder, '/', dataFiles(i).name]);
        
        % Format Data
        switch field

            case 'p'
                PVdata.(plane).positionGrid = content(:,[1,2,3]);
                PVdata.(plane).p.mean = content(:,4);

            case 'U'
                PVdata.(plane).positionGrid = content(:,[1,2,3]);
                PVdata.(plane).u.mean = content(:,4);
                PVdata.(plane).v.mean = content(:,5);
                PVdata.(plane).w.mean = content(:,6);

        end
        
        % Remove Unwanted Planar Deviations
        switch orientation
            
            case 'YZ'
                [count, value] = groupcounts(PVdata.(plane).positionGrid(:,1));
                [~, index] = max(count);
                
                if height(index) == 1
                    index = find(PVdata.(plane).positionGrid(:,1) == value(index));
                else
                    error(['''', plane, ''' Has an Ambiguous Position']);
                end
                
            case 'XZ'
                [count, value] = groupcounts(PVdata.(plane).positionGrid(:,2));
                [~, index] = max(count);
                
                if height(index) == 1
                    index = find(PVdata.(plane).positionGrid(:,2) == value(index));
                else
                    error(['''', plane, ''' Has an Ambiguous Position']);
                end
                
            case 'XY'
                [count, value] = groupcounts(PVdata.(plane).positionGrid(:,3));
                [~, index] = max(count);
                
                if height(index) == 1
                    index = find(PVdata.(plane).positionGrid(:,3) == value(index));
                else
                    error(['''', plane, ''' Has an Ambiguous Position']);
                end
                
        end
        
        PVdata.(plane).positionGrid = PVdata.(plane).positionGrid(index,:);
        
        switch field
            
            case 'p'
                PVdata.(plane).p.mean = PVdata.(plane).p.mean(index);

            case 'U'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
        % Remove Duplicate Entries
        [PVdata.(plane).positionGrid, index] = unique(PVdata.(plane).positionGrid, 'rows', 'stable');
        
        switch field
            
            case 'p'
                PVdata.(plane).p.mean = PVdata.(plane).p.mean(index);

            case 'U'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
        % Sort Position Grid for 'ndgrid' Compatibility
        [PVdata.(plane).positionGrid, index] = sortrows(PVdata.(plane).positionGrid, [3, 2, 1]);
        
        switch field
            
            case 'p'
                PVdata.(plane).p.mean = PVdata.(plane).p.mean(index);

            case 'U'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
    end
    clear i;
    
    PVdata = orderfields(PVdata);

end