%% Planar Probe Data Extractor v1.0
% ----
% Extracts Planar Probe Data From Volumetric Data Processed Using ‘readProbeData.m'
% ----
% Usage: planeData = extractPlanarVelocityProbeData(caseName, volumeData);
%        'caseName'   -> Case Name, Stored as a String
%        'volumeData' -> Volumetric Probe Data, Collated Using 'readProbeData.m'


%% Changelog

% v1.0 - Initial Commit


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Velocity: 'U'


%% Main Function

function planeData = extractPlanarVelocityProbeData(caseName, volumeData)

    disp('Planar Probe Data Extraction');
    disp('-----------------------------');

    disp(' ');
    
    % Temporarily Adjust Data Origin
    if contains(caseName, ["Run_Test", "Windsor"]) && contains(caseName, 'Upstream')
        volumeData.position(:,1) = volumeData.position(:,1) + 1.325;
    end
    
    % Identify Data Limits
    xLimsData = [min(volumeData.position(:,1)); max(volumeData.position(:,1))];
    yLimsData = [min(volumeData.position(:,2)); max(volumeData.position(:,2))];
    zLimsData = [min(volumeData.position(:,3)); max(volumeData.position(:,3))];
    
    disp('Probe Volume:');
    disp(['    X: ', num2str(xLimsData(1)), ' [m] -> ', num2str(xLimsData(2)), ' [m]']);
    disp(['    Y: ', num2str(yLimsData(1)), ' [m] -> ', num2str(yLimsData(2)), ' [m]']);
    disp(['    Z: ', num2str(zLimsData(1)), ' [m] -> ', num2str(zLimsData(2)), ' [m]']);

    disp(' ');

    disp('Possible Planar Orientations:');
    disp('    X (YZ-Plane)');
    disp('    Y (XZ-Plane)');
    disp('    Z (XY-Plane)');
    
    valid = false;
    while ~valid
        disp(' ');
        nPlanes = str2num(input('Select Desired Number of Planes [Scalar]: ', 's')); %#ok<ST2NM>

        if isempty(nPlanes) || isnan(nPlanes) || length(nPlanes) > 1
            disp('    WARNING: Invalid Entry');
        else
            valid = true;
        end

    end
    clear valid;
    
    for i = 1:nPlanes
        
        valid = false;
        while ~valid
            disp(' ');
            selection = input(['Select Orientation of Plane ', num2str(i), ' [X/Y/Z]: '], 's');

            if selection == 'x' | selection == 'X' %#ok<OR2>
                planeData.(['p', num2str(i)]).planeOrientation = 'YZ';
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                planeData.(['p', num2str(i)]).planeOrientation = 'XZ';
                
                valid = true;
            elseif selection == 'z' | selection == 'Z' %#ok<OR2>
                planeData.(['p', num2str(i)]).planeOrientation = 'XY';
                
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;
        
        planeData.(['p', num2str(i)]).xLims = zeros(2,1);
        planeData.(['p', num2str(i)]).yLims = planeData.(['p', num2str(i)]).xLims;
        planeData.(['p', num2str(i)]).zLims = planeData.(['p', num2str(i)]).xLims;

        valid = false;
        while ~valid
            
            switch planeData.(['p', num2str(i)]).planeOrientation
                
                case 'YZ'
                    planeData.(['p', num2str(i)]).yLims = zeros(2,1);
                    planeData.(['p', num2str(i)]).zLims = zeros(2,1);
                    
                    planeData.(['p', num2str(i)]).xLims = inputPos('Planar', 'X');
                    
                    if planeData.(['p', num2str(i)]).xLims == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    
                    if planeData.(['p', num2str(i)]).yLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    
                    if planeData.(['p', num2str(i)]).yLims(2) == -1
                        continue
                    end         
                    
                    planeData.(['p', num2str(i)]).yLims = sort(planeData.(['p', num2str(i)]).yLims);
                    
                    planeData.(['p', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    
                    if planeData.(['p', num2str(i)]).zLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    
                    if planeData.(['p', num2str(i)]).zLims(2) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).zLims = sort(planeData.(['p', num2str(i)]).zLims);
                
                case 'XZ'
                    planeData.(['p', num2str(i)]).xLims = zeros(2,1);
                    planeData.(['p', num2str(i)]).zLims = zeros(2,1);
                    
                    planeData.(['p', num2str(i)]).yLims = inputPos('Planar', 'Y');
                    
                    if planeData.(['p', num2str(i)]).yLims == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    
                    if planeData.(['p', num2str(i)]).xLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    
                    if planeData.(['p', num2str(i)]).xLims(2) == -1
                        continue
                    end         
                    
                    planeData.(['p', num2str(i)]).xLims = sort(planeData.(['p', num2str(i)]).xLims);
                    
                    planeData.(['p', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    
                    if planeData.(['p', num2str(i)]).zLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    
                    if planeData.(['p', num2str(i)]).zLims(2) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).zLims = sort(planeData.(['p', num2str(i)]).zLims);
                    
                case 'XY'
                    planeData.(['p', num2str(i)]).xLims = zeros(2,1);
                    planeData.(['p', num2str(i)]).yLims = zeros(2,1);
                    
                    planeData.(['p', num2str(i)]).zLims = inputPos('Planar', 'Z');
                    
                    if planeData.(['p', num2str(i)]).zLims == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    
                    if planeData.(['p', num2str(i)]).xLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    
                    if planeData.(['p', num2str(i)]).xLims(2) == -1
                        continue
                    end         
                    
                    planeData.(['p', num2str(i)]).xLims = sort(planeData.(['p', num2str(i)]).xLims);
                    
                    planeData.(['p', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    
                    if planeData.(['p', num2str(i)]).yLims(1) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    
                    if planeData.(['p', num2str(i)]).yLims(2) == -1
                        continue
                    end
                    
                    planeData.(['p', num2str(i)]).yLims = sort(planeData.(['p', num2str(i)]).yLims);
            
            end
            
            if (min(planeData.(['p', num2str(i)]).xLims) < xLimsData(1)) || (max(planeData.(['p', num2str(i)]).xLims) > xLimsData(2)) || ...
               (min(planeData.(['p', num2str(i)]).yLims) < yLimsData(1)) || (max(planeData.(['p', num2str(i)]).yLims) > yLimsData(2)) || ...
               (min(planeData.(['p', num2str(i)]).zLims) < zLimsData(1)) || (max(planeData.(['p', num2str(i)]).zLims) > zLimsData(2))
                disp('        WARNING: Plane Lies Outside Probe Volume');
                disp(' ');
            else
                valid = true;
            end
            
        end
        clear valid;
        
    end
    
    disp(' ');
    
    disp('Extracting Planar Probe Data...');
    
    % Shift Requested Plane(s) to Nearest Probe Plane
    for i = 1:nPlanes
        
        switch planeData.(['p', num2str(i)]).planeOrientation
            
            case 'YZ'
                allPlanes = unique(volumeData.position(:,1), 'stable');
                [offset, index] = min(abs(allPlanes - planeData.(['p', num2str(i)]).xLims));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting X: ', num2str(planeData.(['p', num2str(i)]).xLims), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    planeData.(['p', num2str(i)]).xLims = allPlanes(index);
                end
                
                planeData.(['p', num2str(i)]).planePosition = planeData.(['p', num2str(i)]).xLims;
                
                index = find((volumeData.position(:,1) == planeData.(['p', num2str(i)]).xLims) & ...
                             (volumeData.position(:,2) >= planeData.(['p', num2str(i)]).yLims(1) & volumeData.position(:,2) <= planeData.(['p', num2str(i)]).yLims(2)) & ...
                             (volumeData.position(:,3) >= planeData.(['p', num2str(i)]).zLims(1) & volumeData.position(:,3) <= planeData.(['p', num2str(i)]).zLims(2)));
            
            case 'XZ'
                allPlanes = unique(volumeData.position(:,2), 'stable');
                [offset, index] = min(abs(allPlanes - planeData.(['p', num2str(i)]).yLims));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting Y: ', num2str(planeData.(['p', num2str(i)]).yLims), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    planeData.(['p', num2str(i)]).yLims = allPlanes(index);
                end
                
                planeData.(['p', num2str(i)]).planePosition = planeData.(['p', num2str(i)]).yLims;
                
                index = find((volumeData.position(:,1) >= planeData.(['p', num2str(i)]).xLims(1) & volumeData.position(:,1) <= planeData.(['p', num2str(i)]).xLims(2)) & ...
                             (volumeData.position(:,2) == planeData.(['p', num2str(i)]).yLims) & ...
                             (volumeData.position(:,3) >= planeData.(['p', num2str(i)]).zLims(1) & volumeData.position(:,3) <= planeData.(['p', num2str(i)]).zLims(2)));
                
            case 'XY'
                allPlanes = unique(volumeData.position(:,3), 'stable');
                [offset, index] = min(abs(allPlanes - planeData.(['p', num2str(i)]).zLims));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting Z: ', num2str(planeData.(['p', num2str(i)]).zLims), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    planeData.(['p', num2str(i)]).zLims = allPlanes(index);
                end
                
                planeData.(['p', num2str(i)]).planePosition = planeData.(['p', num2str(i)]).zLims;
                
                index = find((volumeData.position(:,1) >= planeData.(['p', num2str(i)]).xLims(1) & volumeData.position(:,1) <= planeData.(['p', num2str(i)]).xLims(2)) & ...
                             (volumeData.position(:,2) >= planeData.(['p', num2str(i)]).yLims(1) & volumeData.position(:,2) <= planeData.(['p', num2str(i)]).yLims(2)) & ...
                             (volumeData.position(:,3) == planeData.(['p', num2str(i)]).zLims));
        
        end
        
        % Collate Planar Probe Data
        planeData.(['p', num2str(i)]).time = volumeData.time;
        planeData.(['p', num2str(i)]).position = volumeData.position(index,:);
        planeData.(['p', num2str(i)]).uMean = volumeData.uMean(index);
        planeData.(['p', num2str(i)]).vMean = volumeData.vMean(index);
        planeData.(['p', num2str(i)]).wMean = volumeData.wMean(index);
        
        planeData.(['p', num2str(i)]).u = cell(height(volumeData.time),1);
        planeData.(['p', num2str(i)]).v =  planeData.(['p', num2str(i)]).u;
        planeData.(['p', num2str(i)]).w =  planeData.(['p', num2str(i)]).u;
        planeData.(['p', num2str(i)]).uPrime = planeData.(['p', num2str(i)]).u;
        planeData.(['p', num2str(i)]).vPrime =  planeData.(['p', num2str(i)]).u;
        planeData.(['p', num2str(i)]).wPrime =  planeData.(['p', num2str(i)]).u;
        
        for j = 1:height(volumeData.time)
            planeData.(['p', num2str(i)]).u{j} = volumeData.u{j}(index,:);
            planeData.(['p', num2str(i)]).v{j} = volumeData.v{j}(index,:);
            planeData.(['p', num2str(i)]).w{j} = volumeData.w{j}(index,:);
            planeData.(['p', num2str(i)]).uPrime{j} = volumeData.uPrime{j}(index,:);
            planeData.(['p', num2str(i)]).vPrime{j} = volumeData.vPrime{j}(index,:);
            planeData.(['p', num2str(i)]).wPrime{j} = volumeData.wPrime{j}(index,:);
        end
        
        % Adjust Limits to Adhere to Data Points
        switch planeData.(['p', num2str(i)]).planeOrientation
            
            case 'YZ'
                planeData.(['p', num2str(i)]).yLims = [min(planeData.(['p', num2str(i)]).position(:,2)); max(planeData.(['p', num2str(i)]).position(:,2))];
                planeData.(['p', num2str(i)]).zLims = [min(planeData.(['p', num2str(i)]).position(:,3)); max(planeData.(['p', num2str(i)]).position(:,3))];
            
            case 'XZ'
                planeData.(['p', num2str(i)]).xLims = [min(planeData.(['p', num2str(i)]).position(:,1)); max(planeData.(['p', num2str(i)]).position(:,1))];
                planeData.(['p', num2str(i)]).zLims = [min(planeData.(['p', num2str(i)]).position(:,3)); max(planeData.(['p', num2str(i)]).position(:,3))];
                
            case 'XY'
                planeData.(['p', num2str(i)]).xLims = [min(planeData.(['p', num2str(i)]).position(:,1)); max(planeData.(['p', num2str(i)]).position(:,1))];
                planeData.(['p', num2str(i)]).yLims = [min(planeData.(['p', num2str(i)]).position(:,2)); max(planeData.(['p', num2str(i)]).position(:,2))];
                
        end
        
        % Revert Data Origin
        if contains(caseName, ["Run_Test", "Windsor"]) && contains(caseName, 'Upstream')
            planeData.(['p', num2str(i)]).position(:,1) = planeData.(['p', num2str(i)]).position(:,1) - 1.325;

            switch planeData.(['p', num2str(i)]).planeOrientation
                
                case 'YZ'
                    planeData.(['p', num2str(i)]).planePosition = planeData.(['p', num2str(i)]).planePosition - 1.325;
                    planeData.(['p', num2str(i)]).xLims = planeData.(['p', num2str(i)]).xLims - 1.325;
                    
                case {'XZ', 'XY'}
                    planeData.(['p', num2str(i)]).xLims = planeData.(['p', num2str(i)]).xLims - 1.325;
            
            end
            
        end
        
        % Rename Plane
        pName = [planeData.(['p', num2str(i)]).planeOrientation, erase(num2str(planeData.(['p', num2str(i)]).planePosition, '%.5f'), '.')];
        planeData = renameStructField(planeData, ['p', num2str(i)], pName);
        planeData.(pName) = orderfields(planeData.(pName));
        
    end
    
    planeData = orderfields(planeData);
    
end


%% Local Functions

function pos = inputPos(type, plane)

    pos = str2double(input(['    ', type, ' ', plane, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        pos = -1;
    end
    
end