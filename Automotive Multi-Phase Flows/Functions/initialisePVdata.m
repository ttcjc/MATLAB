%% ParaView Data Initialisation v2.0
% ----
% Initialisation of Exported ParaView Planar Field Data for Further Processing
% ----
% Usage: [campaignID, caseID, PVdata] = initialisePVdata(saveLoc, field)
%
%        'saveLoc' -> Start of File Path, Stored as a String
%        'field'   -> Desired Field Stored as String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v2.0 - Update To Improve Consistency of Structures Across Repository
% v2.1 - Update To Support Resolved Turbulence Kinetic Energy
% v2.2 - Update To Support Total Pressure Coefficient


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022
% Windsor_Upstream_2023
% Windsor_fullscale


%% Supported Fields

% Velocity: 'U'
% Pressure: 'p'
% Total Pressure Coefficient: 'CpT'
% Resolved TKE: 'kResolved'


%% Main Function

function [campaignID, caseID, PVdata] = initialisePVdata(saveLoc, field)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir([saveLoc, '/Numerical/ParaView'], 'Select Case');
    
    % Confirm Support
    if ~contains(caseFolder, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    campaignID = caseFolder((strfind(caseFolder, 'ParaView/') + 9):(max(strfind(caseFolder, '/')) - 1));
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', caseID]);

    % Confirm Data Availability
    if ~any(strcmp(field, {'U', 'p', 'CpT', 'kResolved'}))
        error('Invalid Field (Available Options: ''p'', ''kResolved'', or ''U'')');
    end

    switch field

        case {'p', 'CpT', 'kResolved'}
            fieldType = 'scalar';
            
        case 'U'
            fieldType = 'vector';

    end
    
    dataFiles = dir([caseFolder, '/', field, '_*.csv']);
    
    if isempty(dataFiles)
        error(['Invalid Case Directory (No ', field, ' Data Found)']);
    end
            
    
    disp(' ');
    disp(' ');

    % Acquire ParaView Data
    disp('ParaView Data Acquisition');
    disp('--------------------------');
    
    disp(' ');

    disp('Initialising:');

    for i = 1:height(dataFiles)
        plane = dataFiles(i).name(1:(end - 4)); % Ignore '.csv'

        disp(['    ', plane]);
        
        % Confirm Valid Data Orientation
        if contains(plane, 'Base') || contains(plane, '_X_')
            orientation = 'YZ';
        elseif contains(plane, '_Y_')
            orientation = 'XZ';
        elseif contains(plane, '_Z_')
            orientation = 'XY';
        else
            error(['''', plane, ''' Is an Invalid Data File (Unexpected Naming Convention)']);
        end
    
        % Load Data
        content = importdata([caseFolder, '/', dataFiles(i).name]);
        
        % Sort Columns
        index = find(strcmp(content.colheaders, '"Points:0"'));
        
        if index ~= 1
            
            switch fieldType
                
                case 'scalar'
                    content = content.data(:,[2,3,4,1]);
                    
                case 'vector'
                    content = content.data(:,[4,5,6,1,2,3]);
                    
            end
            
        else
            content = content.data;
        end
        
        % Format Data
        switch fieldType

            case 'scalar'
                PVdata.(plane).positionGrid = content(:,[1,2,3]);
                PVdata.(plane).(field).mean = content(:,4);

            case 'vector'
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
        
        switch fieldType
            
            case 'scalar'
                PVdata.(plane).(field).mean = PVdata.(plane).(field).mean(index);

            case 'vector'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
        % Remove Duplicate Entries
        [PVdata.(plane).positionGrid, index] = unique(PVdata.(plane).positionGrid, 'rows', 'stable');
        
        switch fieldType
            
            case 'scalar'
                PVdata.(plane).(field).mean = PVdata.(plane).(field).mean(index);

            case 'vector'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
        % Sort Position Grid for 'ndgrid' Compatibility
        [PVdata.(plane).positionGrid, index] = sortrows(PVdata.(plane).positionGrid, [3, 2, 1]);
        
        switch fieldType
            
            case 'scalar'
                PVdata.(plane).(field).mean = PVdata.(plane).(field).mean(index);

            case 'vector'
                PVdata.(plane).u.mean = PVdata.(plane).u.mean(index);
                PVdata.(plane).v.mean = PVdata.(plane).v.mean(index);
                PVdata.(plane).w.mean = PVdata.(plane).w.mean(index);
        
        end
        
    end
    clear i;
    
    PVdata = orderfields(PVdata);

end
