%% Experimental Flow Data Initialisation v2.0
% ----
% Initialisation of Experimental Flow Data for Further Processing
% ----
% Usage: [campaignID, caseID, expData] = initialiseExpFlowData(saveLocation, field);
%
%        'saveLocation' -> Start of File Path, Stored as a String
%        'field'        -> Desired Field Stored as String
%        'normDims'     -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Update To Improve Consistency of Structures Across Repository


%% Supported Experimental Cases

% Windsor Model Experimental Campaign (Varney)


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [campaignID, caseID, expData] = initialiseExpFlowData(saveLocation, field)

    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir([saveLocation, '/Experimental'], 'Select Case');
    
    % Confirm Support
    if ~contains(caseFolder, 'Windsor Model (Varney)')
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    campaignID = 'Varney';
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', caseFolder]);


    % Confirm Data Availability
    if ~strcmp(field, 'p') && ~strcmp(field, 'U')
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end
    
    switch field

        case 'p'
            fieldLabel = 'Pressure';

            testFiles = dir([caseFolder, '/', fieldLabel, '/*.mat']);

        case 'U'
            fieldLabel = 'PIV';

            testFiles = dir([caseFolder, '/', fieldLabel, '/*.mat']);

    end
    
    if isempty(testFiles)
        error(['Invalid Case Directory (No ', fieldLabel, ' Data Found)']);
    end
    
    disp(' ');
    disp(' ');
    
    % Acquire ParaView Data
    disp('Experimental Data Acquisition');
    disp('------------------------------');

    disp(' ');
    
    disp('Initialising:');

    for i = 1:height(testFiles)
        testID = testFiles(i).name(1:(end - 4)); % Ignore '.mat'
        
        disp(['    ', testID]);
        
        % Load Data
        content = load([caseFolder, '/', fieldLabel, '/', testFiles(i).name]);

        % Format Data                
        switch field

            case 'p'
                expData.(testID).positionGrid = content.xyz(4).values(:,[1,2,3]) / 1000; % [mm -> m]
                expData.(testID).positionGrid(:,3) = -expData.(testID).positionGrid(:,3);

                expData.(testID).Cp.mean = content.Cpmean(4).values;
                expData.(testID).Cp.RMS = content.Cprms(4).values;

            case 'U'
                expData.(testID).positionGrid(:,1) = content.x.values(:) / 1000;
                expData.(testID).positionGrid(:,2) = content.y.values(:) / 1000;
                expData.(testID).positionGrid(:,3) = content.z.values(:) / 1000;

                if contains(testID, '_X_')
                    expData.(testID).positionGrid(:,1) = str2double([testID(end - 6), '.', testID((end - 4):end)]);
                    
                    expData.(testID).u.mean = double(content.Vx_mean.values(:));
                    expData.(testID).v.mean = double(content.Vy_mean.values(:));
                    expData.(testID).w.mean = double(content.Vz_mean.values(:));
                    
                    expData.(testID).u.RMS = double(content.Vx_std.values(:));
                    expData.(testID).v.RMS = double(content.Vy_std.values(:));
                    expData.(testID).w.RMS = double(content.Vz_std.values(:));
                elseif contains(testID, '_Y_')
                    expData.(testID).positionGrid(:,1) = -expData.(testID).positionGrid(:,1);
                    expData.(testID).positionGrid(:,2) = str2double([testID(end - 6), '.', testID((end - 4):end)]);
                    
                    expData.(testID).u.mean = double(-content.Vx_mean.values(:));
                    expData.(testID).v.mean = nan(height(expData.(testID).positionGrid),1);
                    expData.(testID).w.mean = double(content.Vz_mean.values(:));
                    
                    expData.(testID).u.RMS = double(content.Vx_std.values(:));
                    expData.(testID).v.RMS = nan(height(expData.(testID).positionGrid),1);
                    expData.(testID).w.RMS = double(content.Vz_std.values(:));
                elseif contains(testID, '_Z_')
                    expData.(testID).positionGrid(:,3) = str2double([testID(end - 6), '.', testID((end - 4):end)]);
                    
                    expData.(testID).u.mean = double(content.Vx_mean.values(:));
                    expData.(testID).v.mean = double(content.Vy_mean.values(:));
                    expData.(testID).w.mean = nan(height(expData.(testID).positionGrid),1);
                    
                    expData.(testID).u.RMS = double(content.Vx_std.values(:));
                    expData.(testID).v.RMS = double(content.Vy_std.values(:));
                    expData.(testID).w.RMS = nan(height(expData.(testID).positionGrid),1);
                else
                    error('Invalid Test File (No Orientation Information Available)');
                end

        end
        
        % Remove NaN Values
        switch field

            case 'p'
                index = find(isnan(expData.(testID).Cp.mean));

                expData.(testID).positionGrid(index,:) = [];
                expData.(testID).Cp.mean(index) = [];
                expData.(testID).Cp.RMS(index) = [];                        

            case 'U'

                if contains(testID, '_X_')
                    index = find(isnan(expData.(testID).u.mean) | isnan(expData.(testID).v.mean) | ...
                                 isnan(expData.(testID).w.mean));
                elseif contains(testID, '_Y_')
                    index = find(isnan(expData.(testID).u.mean) | isnan(expData.(testID).w.mean));
                elseif contains(testID, '_Z_')
                    index = find(isnan(expData.(testID).u.mean) | isnan(expData.(testID).v.mean));
                end

                expData.(testID).positionGrid(index,:) = [];
                expData.(testID).u.mean(index) = [];
                expData.(testID).v.mean(index) = [];
                expData.(testID).w.mean(index) = [];
                expData.(testID).u.RMS(index) = [];
                expData.(testID).v.RMS(index) = [];
                expData.(testID).w.RMS(index) = [];

        end
        
        % Remove Duplicate Entries
        [expData.(testID).positionGrid, index] = unique(expData.(testID).positionGrid, 'stable', 'rows');

        switch field

            case 'p'
                expData.(testID).Cp.mean = expData.(testID).Cp.mean(index);
                expData.(testID).Cp.RMS = expData.(testID).Cp.RMS(index);

            case 'U'
                expData.(testID).u.mean = expData.(testID).u.mean(index);
                expData.(testID).v.mean = expData.(testID).v.mean(index);
                expData.(testID).w.mean = expData.(testID).w.mean(index);
                expData.(testID).u.RMS = expData.(testID).u.RMS(index);
                expData.(testID).v.RMS = expData.(testID).v.RMS(index);
                expData.(testID).w.RMS = expData.(testID).w.RMS(index);

        end
        
        % Adjust Data Origin
        switch field

            case 'p'
                expData.(testID).positionGrid(:,1) = expData.(testID).positionGrid(:,1) + ...
                                                 (0.48325 - min(expData.(testID).positionGrid(:,1)));

            case 'U'

                if contains(testID, '_X_')
                    % Adjust
                elseif contains(testID, '_Y_')

                    if contains(testID, '040341')
                        expData.(testID).positionGrid(:,1) = expData.(testID).positionGrid(:,1) + ...
                                                             (0.48325 - min(expData.(testID).positionGrid(:,1))) - 0.032;
                        expData.(testID).positionGrid(:,3) = expData.(testID).positionGrid(:,3) + ...
                                                             (0.029 - min(expData.(testID).positionGrid(:,3))) - 0.0265;
                    else
                        expData.(testID).positionGrid(:,1) = expData.(testID).positionGrid(:,1) + ...
                                                             (0.48325 - min(expData.(testID).positionGrid(:,1))) - 0.004;
                        expData.(testID).positionGrid(:,3) = expData.(testID).positionGrid(:,3) + ...
                                                             (0.029 - min(expData.(testID).positionGrid(:,3)));
                    end

                elseif contains(testID, '_Z_')
                    expData.(testID).positionGrid(:,1) = expData.(testID).positionGrid(:,1) + ...
                                                         (0.48325 - min(expData.(testID).positionGrid(:,1))) - 0.008;
                end

        end
        
        % Sort Position Grid for 'ndgrid' Compatibility
        [expData.(testID).positionGrid, index] = sortrows(expData.(testID).positionGrid, [3, 2, 1]);

        switch field

            case 'p'
                expData.(testID).Cp.mean = expData.(testID).Cp.mean(index);
                expData.(testID).Cp.RMS = expData.(testID).Cp.RMS(index);

            case 'U'
                expData.(testID).u.mean = expData.(testID).u.mean(index);
                expData.(testID).v.mean = expData.(testID).v.mean(index);
                expData.(testID).w.mean = expData.(testID).w.mean(index);
                expData.(testID).u.RMS = expData.(testID).u.RMS(index);
                expData.(testID).v.RMS = expData.(testID).v.RMS(index);
                expData.(testID).w.RMS = expData.(testID).w.RMS(index);

        end
        
    end
    clear i;
    
    expData = orderfields(expData);
    
end