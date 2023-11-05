%% Experimental Data Initialisation v1.0
% ----
% Initialisation of OpenFOAM v7 Probe Data for Further Processing
% ----
% Usage: [campaign, expData, geometry, ...
%         xDims, yDims, zDims, spacePrecision] = initialiseExpData(field, normalise);
%
%        'field'     -> Desired Field Stored as String
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit


%% Supported Experimental Cases

% Windsor Model Experimental Campaign (Varney)


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [caseName, expData, geometry, ...
          xDims, yDims, zDims, spacePrecision] = initialiseExpData(field, normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir('~/Data', 'Select Case');

    disp(' ');

    disp(['Case: ', caseFolder]);

    % Confirm Support
    if contains(caseFolder, 'Windsor Model (Varney)')
        caseName = 'Varney';
    else
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    % Confirm Data Availability
    if ~strcmp(field, 'p') && ~strcmp(field, 'U')
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end

    switch caseName

        case 'Varney'

            switch field
        
                case 'p'
                    fieldLabel = 'Pressure';

                    if exist([caseFolder, '/', fieldLabel], 'dir')
                        testFiles = dir([caseFolder, '/', fieldLabel, '/*.mat']);
                    else
                        error(['Invalid Case Directory (No ', fieldLabel, ' Data Found)']);
                    end
                
                case 'U'
                    fieldLabel = 'PIV';

                    if exist([caseFolder, '/', fieldLabel], 'dir')
                        testFiles = dir([caseFolder, '/', fieldLabel, '/*.mat']);
                    else
                        error(['Invalid Case Directory (No ', fieldLabel, ' Data Found)']);
                    end
        
            end
        
        if isempty(testFiles)
            error(['Invalid Case Directory (No ', fieldLabel, ' Data Found)']);
        end

    end
    
    disp(' ');
    disp(' ');
    
    % Acquire ParaView Data
    disp('Experimental Data Acquisition');
    disp('------------------------------');

    disp(' ');
    
    disp('Formatting:');

    for i = 1:height(testFiles)
        testID = testFiles(i).name(1:(end - 4)); % Ignore .mat
        disp(['    ', testID]);
        
        % Load Data
        content = load([caseFolder, '/', fieldLabel, '/', testFiles(i).name]);

        % Format Data
        switch caseName

            case 'Varney'
                
                switch field
        
                    case 'p'
                        expData.(testID).position = content.xyz(4).values(:,[1,2,3]) / 1000; % [mm -> m]
                        expData.(testID).position(:,3) = -expData.(testID).position(:,3);
                        
                        expData.(testID).CpMean = content.Cpmean(4).values;
                        expData.(testID).CpRMS = content.Cprms(4).values;

                        expData.(testID).planeOrientation = 'YZ';
                    
                    case 'U'
                        expData.(testID).position(:,1) = content.x.values(:) / 1000;
                        expData.(testID).position(:,2) = content.y.values(:) / 1000;
                        expData.(testID).position(:,3) = content.z.values(:) / 1000;
                        
                        if contains(testID, '_x')
                            expData.(testID).uMean = double(content.Vx_mean.values(:));
                            expData.(testID).vMean = double(content.Vy_mean.values(:));
                            expData.(testID).wMean = double(content.Vz_mean.values(:));
                            expData.(testID).uRMS = double(content.Vx_std.values(:));
                            expData.(testID).vRMS = double(content.Vy_std.values(:));
                            expData.(testID).wRMS = double(content.Vz_std.values(:));

                            expData.(testID).planeOrientation = 'YZ';
                            expData.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            expData.(testID).position(:,1) = expData.(testID).planePosition;
                        elseif contains(testID, '_y')
                            expData.(testID).uMean = double(-content.Vx_mean.values(:));
                            expData.(testID).vMean = nan(height(expData.(testID).position),1);
                            expData.(testID).wMean = double(content.Vz_mean.values(:));
                            expData.(testID).uRMS = double(content.Vx_std.values(:));
                            expData.(testID).vRMS = nan(height(expData.(testID).position),1);
                            expData.(testID).wRMS = double(content.Vz_std.values(:));

                            expData.(testID).planeOrientation = 'XZ';
                            expData.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            expData.(testID).position(:,1) = -expData.(testID).position(:,1);
                            expData.(testID).position(:,2) = expData.(testID).planePosition;
                        elseif contains(testID, '_z')
                            expData.(testID).uMean = double(content.Vx_mean.values(:));
                            expData.(testID).vMean = double(content.Vy_mean.values(:));
                            expData.(testID).wMean = nan(height(expData.(testID).position),1);
                            expData.(testID).uRMS = double(content.Vx_std.values(:));
                            expData.(testID).vRMS = double(content.Vy_std.values(:));
                            expData.(testID).wRMS = nan(height(expData.(testID).position),1);

                            expData.(testID).planeOrientation = 'XY';
                            expData.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            expData.(testID).position(:,3) = expData.(testID).planePosition;
                        else
                            error('Invalid Test File (No Orientation Information Available)');
                        end
                
                end
        
        end

        % Remove NaN Values
        switch caseName

            case 'Varney'

                switch field

                    case 'p'
                        index = find(isnan(expData.(testID).CpMean));

                        expData.(testID).position(index,:) = [];
                        expData.(testID).CpMean(index) = [];
                        expData.(testID).CpRMS(index) = [];                        

                    case 'U'

                        if contains(testID, '_x')
                            index = find(isnan(expData.(testID).vMean));
                        elseif contains(testID, '_y')
                            index = find(isnan(expData.(testID).wMean));
                        elseif contains(testID, '_z')
                            index = find(isnan(expData.(testID).uMean));
                        end

                        expData.(testID).position(index,:) = [];
                        expData.(testID).uMean(index) = [];
                        expData.(testID).vMean(index) = [];
                        expData.(testID).wMean(index) = [];
                        expData.(testID).uRMS(index) = [];
                        expData.(testID).vRMS(index) = [];
                        expData.(testID).wRMS(index) = [];
                
                end

        end

        % Remove Duplicate Entries
        switch caseName

            case 'Varney'
                [expData.(testID).position, index] = unique(expData.(testID).position, 'stable', 'rows');

                switch field

                    case 'p'
                        expData.(testID).CpMean = expData.(testID).CpMean(index);
                        expData.(testID).CpRMS = expData.(testID).CpRMS(index);

                    case 'U'
                        expData.(testID).uMean = expData.(testID).uMean(index);
                        expData.(testID).vMean = expData.(testID).vMean(index);
                        expData.(testID).wMean = expData.(testID).wMean(index);
                        expData.(testID).uRMS = expData.(testID).uRMS(index);
                        expData.(testID).vRMS = expData.(testID).vRMS(index);
                        expData.(testID).wRMS = expData.(testID).wRMS(index);

                end

        end
        
        % Adjust Data Origin
        switch caseName

            case 'Varney'
                
                switch field

                    case 'p'
                        expData.(testID).position(:,1) = expData.(testID).position(:,1) + ...
                                                         (0.48325 - min(expData.(testID).position(:,1)));
                                                  
                        expData.(testID).planePosition = expData.(testID).position(1,1);

                    case 'U'
                        
                        if contains(testID, '_x')
                            % Adjust
                        elseif contains(testID, '_y')
                            
                            if contains(testID, '040341')
                                expData.(testID).position(:,1) = expData.(testID).position(:,1) + ...
                                                                 (0.48325 - min(expData.(testID).position(:,1))) - 0.032;
                                expData.(testID).position(:,3) = expData.(testID).position(:,3) + ...
                                                                 (0.029 - min(expData.(testID).position(:,3))) - 0.0265;
                            else
                                expData.(testID).position(:,1) = expData.(testID).position(:,1) + ...
                                                                 (0.48325 - min(expData.(testID).position(:,1))) - 0.004;
                                expData.(testID).position(:,3) = expData.(testID).position(:,3) + ...
                                                                 (0.029 - min(expData.(testID).position(:,3)));
                            end
                                                  
                        elseif contains(testID, '_z')
                            expData.(testID).position(:,1) = expData.(testID).position(:,1) + ...
                                                             (0.48325 - min(expData.(testID).position(:,1))) - 0.008;
                        end
                        
                end

        end
        
        % Crop Data Boundaries
        switch caseName

            case 'Varney'
                
                switch field

                    case 'p'
                        % Crop

                    case 'U'
                        
                        if contains(testID, '_x')
                            % Crop
                        elseif contains(testID, '_y')
                            index = find(expData.(testID).position(:,1) >= 0.48325 & (expData.(testID).position(:,3) >= 0.035 & expData.(testID).position(:,2) <= 0.354));
                            
                            expData.(testID).position = expData.(testID).position(index,:);
                            expData.(testID).uMean = expData.(testID).uMean(index);
                            expData.(testID).vMean = expData.(testID).vMean(index);
                            expData.(testID).wMean = expData.(testID).wMean(index);
                            expData.(testID).uRMS = expData.(testID).uRMS(index);
                            expData.(testID).vRMS = expData.(testID).vRMS(index);
                            expData.(testID).wRMS = expData.(testID).wRMS(index);
                                                  
                        elseif contains(testID, '_z')
                            index = find(expData.(testID).position(:,1) >= 0.48325 & (expData.(testID).position(:,2) >= -0.21 & expData.(testID).position(:,2) <= 0.21));
                            
                            expData.(testID).position = expData.(testID).position(index,:);
                            expData.(testID).uMean = expData.(testID).uMean(index);
                            expData.(testID).vMean = expData.(testID).vMean(index);
                            expData.(testID).wMean = expData.(testID).wMean(index);
                            expData.(testID).uRMS = expData.(testID).uRMS(index);
                            expData.(testID).vRMS = expData.(testID).vRMS(index);
                            expData.(testID).wRMS = expData.(testID).wRMS(index);
                        end
                        
                end

        end
        
    end
    clear i;

    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);
    
    % Normalise Data Dimensions
    for i = 1:height(testFiles)
        testID = testFiles(i).name(1:(end - 4));
        
        if normalise

            switch caseName

            case 'Varney'
                expData.(testID).position = round((expData.(testID).position / 1.044), spacePrecision);
                expData.(testID).planePosition = round((expData.(testID).planePosition / 1.044), spacePrecision);
            end

        end
        
        expData.(testID) = orderfields(expData.(testID));
    end
    clear i;
    
    expData = orderfields(expData);
    
end