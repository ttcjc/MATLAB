%% Experimental Data Initialisation v1.0
% ----
% Initialisation of OpenFOAM v7 Probe Data for Further Processing
% ----
% Usage: [campaign, data, geometry, xDims, yDims, zDims, precision] = initialiseExpData(field, normalise);
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

function [campaign, data, geometry, xDims, yDims, zDims, precision] = initialiseExpData(field, normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir('~/Data/Experimental', 'Select Case');

    disp(' ');

    disp(['Case: ', caseFolder]);

    % Confirm Support
    if contains(caseFolder, 'Windsor Model (Varney)')
        campaign = 'Varney';
    else
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    % Confirm Data Availability
    if ~strcmp(field, 'p') && ~strcmp(field, 'U')
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end

    switch campaign

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
        switch campaign

            case 'Varney'
                
                switch field
        
                    case 'p'
                        data.(testID).position = content.xyz(4).values(:,[1,2,3]) / 1000; % [mm -> m]

                        data.(testID).CpMean = content.Cpmean(4).values;
                        data.(testID).CpRMS = content.Cprms(4).values;

                        data.(testID).planeOrientation = 'X';
                    
                    case 'U'
                        data.(testID).position(:,1) = content.x.values(:) / 1000;
                        data.(testID).position(:,2) = content.y.values(:) / 1000;
                        data.(testID).position(:,3) = content.z.values(:) / 1000;
                        
                        if contains(testID, '_x')
                            data.(testID).uMean = double(content.Vx_mean.values(:));
                            data.(testID).vMean = double(content.Vy_mean.values(:));
                            data.(testID).wMean = double(content.Vz_mean.values(:));
                            data.(testID).uRMS = double(content.Vx_std.values(:));
                            data.(testID).vRMS = double(content.Vy_std.values(:));
                            data.(testID).wRMS = double(content.Vz_std.values(:));

                            data.(testID).planeOrientation = 'X';
                            data.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            data.(testID).position(:,1) = data.(testID).planePosition;
                        elseif contains(testID, '_y')
                            data.(testID).uMean = double(-content.Vx_mean.values(:));
                            data.(testID).vMean = nan(height(data.(testID).position),1);
                            data.(testID).wMean = double(content.Vz_mean.values(:));
                            data.(testID).uRMS = double(content.Vx_std.values(:));
                            data.(testID).vRMS = nan(height(data.(testID).position),1);
                            data.(testID).wRMS = double(content.Vz_std.values(:));

                            data.(testID).planeOrientation = 'Y';
                            data.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            data.(testID).position(:,1) = -data.(testID).position(:,1);
                            data.(testID).position(:,2) = data.(testID).planePosition;
                        elseif contains(testID, '_z')
                            data.(testID).uMean = double(content.Vx_mean.values(:));
                            data.(testID).vMean = double(content.Vy_mean.values(:));
                            data.(testID).wMean = nan(height(data.(testID).position),1);
                            data.(testID).uRMS = double(content.Vx_std.values(:));
                            data.(testID).vRMS = double(content.Vy_std.values(:));
                            data.(testID).wRMS = nan(height(data.(testID).position),1);

                            data.(testID).planeOrientation = 'Z';
                            data.(testID).planePosition = str2double([testID(end - 5), '.', testID((end - 4):end)]);
                            data.(testID).position(:,3) = data.(testID).planePosition;
                        else
                            error('Invalid Test File (No Orientation Information Available)');
                        end
                
                end
        
        end

        % Remove NaN Values
        switch campaign

            case 'Varney'

                switch field

                    case 'p'
                        index = find(isnan(data.(testID).CpMean));

                        data.(testID).position(index,:) = [];
                        data.(testID).CpMean(index) = [];
                        data.(testID).CpRMS(index) = [];                        

                    case 'U'

                        if contains(testID, '_x')
                            index = find(isnan(data.(testID).vMean));
                        elseif contains(testID, '_y')
                            index = find(isnan(data.(testID).wMean));
                        elseif contains(testID, '_z')
                            index = find(isnan(data.(testID).uMean));
                        end

                        data.(testID).position(index,:) = [];
                        data.(testID).uMean(index) = [];
                        data.(testID).vMean(index) = [];
                        data.(testID).wMean(index) = [];
                        data.(testID).uRMS(index) = [];
                        data.(testID).vRMS(index) = [];
                        data.(testID).wRMS(index) = [];
                
                end

        end

        % Remove Duplicate Entries
        switch campaign

            case 'Varney'
                [data.(testID).position, index] = unique(data.(testID).position, 'stable', 'rows');

                switch field

                    case 'p'
                        data.(testID).CpMean = data.(testID).CpMean(index);
                        data.(testID).CpRMS = data.(testID).CpRMS(index);

                    case 'U'
                        data.(testID).uMean = data.(testID).uMean(index);
                        data.(testID).vMean = data.(testID).vMean(index);
                        data.(testID).wMean = data.(testID).wMean(index);
                        data.(testID).uRMS = data.(testID).uRMS(index);
                        data.(testID).vRMS = data.(testID).vRMS(index);
                        data.(testID).wRMS = data.(testID).wRMS(index);

                end

        end
        
        % Adjust Data Origin
        switch campaign

            case 'Varney'
                
                switch field

                    case 'p'
                        data.(testID).position(:,1) = data.(testID).position(:,1) + ...
                                                      (0.48325 - min(data.(testID).position(:,1)));
                                                  
                        data.(testID).planePosition = data.(testID).position(1,1);

                    case 'U'
                        data.(testID).position(:,1) = data.(testID).position(:,1) + ...
                                                      (0.48325 - min(data.(testID).position(:,1)));
                        data.(testID).position(:,3) = data.(testID).position(:,3) + ...
                                                      (0.029 - min(data.(testID).position(:,3)));

                end

        end
        
    end

    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, precision] = selectGeometry(normalise);
    
    % Normalise Data Dimensions
    for i = 1:height(testFiles)
        testID = testFiles(i).name(1:(end - 4));
        
        if normalise

            switch campaign

            case 'Varney'
                data.(testID).position(:,1) = round(data.(testID).position(:,1) / 1.044, precision);
                data.(testID).position(:,2) = round(data.(testID).position(:,2) / 1.044, precision);
                data.(testID).position(:,3) = round(data.(testID).position(:,3) / 1.044, precision);

                data.(testID).planePosition = round(data.(testID).planePosition / 1.044, precision);

            end

        end
        
        data.(testID) = orderfields(data.(testID));
    end
    
    data = orderfields(data);
    
end
