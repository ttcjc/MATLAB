%% Planar Probe Data Extractor v1.0
% ----
% Extracts Planar Probe Data From Volumetric Data Processed Using ‘readProbeData.m'
% ----
% Usage: data = extractPlanarProbeData();


%% Changelog

% v1.0 - Initial Commit


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Velocity: 'U'


%% Main Function

function data = extractPlanarProbeData(caseFolder, vData)

    disp('Planar Probe Data Extraction');
    disp('-----------------------------');

    disp(' ');
    
    % Temporarily Adjust Data Origin
    if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
        vData.position(:,1) = vData.position(:,1) + 1.325;
    end
    
    % Identify Data Limits
    xLimsData = [min(vData.position(:,1)); max(vData.position(:,1))];
    yLimsData = [min(vData.position(:,2)); max(vData.position(:,2))];
    zLimsData = [min(vData.position(:,3)); max(vData.position(:,3))];
    
    disp('Probe Volume:');
    disp(['    X: ', num2str(xLimsData(1)), ' [m] -> ', num2str(xLimsData(2)), ' [m]']);
    disp(['    Y: ', num2str(yLimsData(1)), ' [m] -> ', num2str(yLimsData(2)), ' [m]']);
    disp(['    Z: ', num2str(zLimsData(1)), ' [m] -> ', num2str(zLimsData(2)), ' [m]']);

    disp(' ');

    disp('Possible Planar Orientations:');
    disp('    X: Normal [1 0 0]');
    disp('    Y: Normal [0 1 0]');
    disp('    Z: Normal [0 0 1]');
    
    valid = false;
    while ~valid
        disp(' ');
        nPlanes = str2num(input('Select Desired Number of Planes [Scalar]: ', 's')); %#ok<ST2NM>

        if isnan(nPlanes) || length(nPlanes) > 1
            disp('    WARNING: Invalid Entry');
        else
            valid = true;
        end

    end
    
    for i = 1:nPlanes
        
        valid = false;
        while ~valid
            disp(' ');
            selection = input(['Select Orientation of Plane ', num2str(i), ' [X/Y/Z]: '], 's');

            if selection == 'x' | selection == 'X' %#ok<OR2>
                data.(['p', num2str(i)]).planeOrientation = 'X';
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                data.(['p', num2str(i)]).planeOrientation = 'Y';
                valid = true;
            elseif selection == 'z' | selection == 'Z' %#ok<OR2>
                data.(['p', num2str(i)]).planeOrientation = 'Z';
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
        data.(['p', num2str(i)]).xLims = zeros(2,1);
        data.(['p', num2str(i)]).yLims = data.(['p', num2str(i)]).xLims;
        data.(['p', num2str(i)]).zLims = data.(['p', num2str(i)]).xLims;

        valid = false;
        while ~valid
            
            switch data.(['p', num2str(i)]).planeOrientation
                
                case 'X'
                    data.(['p', num2str(i)]).xLims(1) = inputPos('Planar', 'X');
                    data.(['p', num2str(i)]).xLims(2) = data.(['p', num2str(i)]).xLims(1);

                    data.(['p', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    data.(['p', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    data.(['p', num2str(i)]).yLims = sort(data.(['p', num2str(i)]).yLims);

                    data.(['p', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    data.(['p', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    data.(['p', num2str(i)]).zLims = sort(data.(['p', num2str(i)]).zLims);
                    
                case 'Y'
                    data.(['p', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    data.(['p', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    data.(['p', num2str(i)]).xLims = sort(data.(['p', num2str(i)]).xLims);

                    data.(['p', num2str(i)]).yLims(1) = inputPos('Planar', 'Y');
                    data.(['p', num2str(i)]).yLims(2) = data.(['p', num2str(i)]).yLims(1);

                    data.(['p', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    data.(['p', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    data.(['p', num2str(i)]).zLims = sort(data.(['p', num2str(i)]).zLims);
                case 'Z'
                    data.(['p', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    data.(['p', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    data.(['p', num2str(i)]).xLims = sort(data.(['p', num2str(i)]).xLims);

                    data.(['p', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    data.(['p', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    data.(['p', num2str(i)]).yLims = sort(data.(['p', num2str(i)]).yLims);

                    data.(['p', num2str(i)]).zLims(1) = inputPos('Planar', 'Z');
                    data.(['p', num2str(i)]).zLims(2) = data.(['p', num2str(i)]).zLims(1);
            
            end
            
            if data.(['p', num2str(i)]).xLims(1) < xLimsData(1) || data.(['p', num2str(i)]).xLims(2) > xLimsData(2) || ...
               data.(['p', num2str(i)]).yLims(1) < yLimsData(1) || data.(['p', num2str(i)]).yLims(2) > yLimsData(2) || ...
               data.(['p', num2str(i)]).zLims(1) < zLimsData(1) || data.(['p', num2str(i)]).zLims(2) > zLimsData(2)
                disp('        WARNING: Plane Lies Outside Probe Volume');
                disp(' ');
            else
                valid = true;
            end
            
        end
        
    end
    
    disp(' ');
    
    disp('Extracting Planar Probe Data...');
    
    % Shift Requested Plane(s) to Nearest Probe Plane
    for i = 1:nPlanes
        
        switch data.(['p', num2str(i)]).planeOrientation
            
            case 'X'
                allPlanes = unique(vData.position(:,1), 'stable');
                [offset, index] = min(abs(allPlanes - data.(['p', num2str(i)]).xLims(1)));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting X: ', num2str(data.(['p', num2str(i)]).xLims(1)), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    data.(['p', num2str(i)]).xLims(1) = allPlanes(index);
                    data.(['p', num2str(i)]).xLims(2) = data.(['p', num2str(i)]).xLims(1);
                end
                
                data.(['p', num2str(i)]).planePosition = data.(['p', num2str(i)]).xLims(1);
                
                index = find((vData.position(:,1) == data.(['p', num2str(i)]).xLims(1)) & ...
                             (vData.position(:,2) >= data.(['p', num2str(i)]).yLims(1) & vData.position(:,2) <= data.(['p', num2str(i)]).yLims(2)) & ...
                             (vData.position(:,3) >= data.(['p', num2str(i)]).zLims(1) & vData.position(:,3) <= data.(['p', num2str(i)]).zLims(2)));
            
            case 'Y'
                allPlanes = unique(vData.position(:,2), 'stable');
                [offset, index] = min(abs(allPlanes - data.(['p', num2str(i)]).yLims(1)));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting Y: ', num2str(data.(['p', num2str(i)]).yLims(1)), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    data.(['p', num2str(i)]).yLims(1) = allPlanes(index);
                    data.(['p', num2str(i)]).yLims(2) = data.(['p', num2str(i)]).yLims(1);
                end
                
                data.(['p', num2str(i)]).planePosition = data.(['p', num2str(i)]).yLims(1);
                
                index = find((vData.position(:,1) >= data.(['p', num2str(i)]).xLims(1) & vData.position(:,1) <= data.(['p', num2str(i)]).xLims(2)) & ...
                             (vData.position(:,2) == data.(['p', num2str(i)]).yLims(2)) & ...
                             (vData.position(:,3) >= data.(['p', num2str(i)]).zLims(1) & vData.position(:,3) <= data.(['p', num2str(i)]).zLims(2)));
                
            case 'Z'
                allPlanes = unique(vData.position(:,3), 'stable');
                [offset, index] = min(abs(allPlanes - data.(['p', num2str(i)]).zLims(1)));
                
                if offset ~= 0
                    disp('    NOTE: Requested Plane Unavailable');
                    disp(['        Shifting Z: ', num2str(data.(['p', num2str(i)]).zLims(1)), ' [m] -> ', num2str(allPlanes(index)), ' [m]']);
                    data.(['p', num2str(i)]).zLims(1) = allPlanes(index);
                    data.(['p', num2str(i)]).zLims(2) = data.(['p', num2str(i)]).zLims(1);
                end
                
                data.(['p', num2str(i)]).planePosition = data.(['p', num2str(i)]).zLims(1);
                
                index = find((vData.position(:,1) >= data.(['p', num2str(i)]).xLims(1) & vData.position(:,1) <= data.(['p', num2str(i)]).xLims(2)) & ...
                             (vData.position(:,2) >= data.(['p', num2str(i)]).yLims(1) & vData.position(:,2) <= data.(['p', num2str(i)]).yLims(2)) & ...
                             (vData.position(:,3) == data.(['p', num2str(i)]).zLims(2)));
        
        end
        
        % Collate Planar Probe Data
        data.(['p', num2str(i)]).time = vData.time;
        data.(['p', num2str(i)]).position = vData.position(index,:);
        data.(['p', num2str(i)]).uMean = vData.uMean(index);
        data.(['p', num2str(i)]).vMean = vData.vMean(index);
        data.(['p', num2str(i)]).wMean = vData.wMean(index);
        
        data.(['p', num2str(i)]).u = cell(height(vData.time),1);
        data.(['p', num2str(i)]).v =  data.(['p', num2str(i)]).u;
        data.(['p', num2str(i)]).w =  data.(['p', num2str(i)]).u;
        data.(['p', num2str(i)]).uPrime = data.(['p', num2str(i)]).u;
        data.(['p', num2str(i)]).vPrime =  data.(['p', num2str(i)]).u;
        data.(['p', num2str(i)]).wPrime =  data.(['p', num2str(i)]).u;
        
        for j = 1:height(vData.time)
            data.(['p', num2str(i)]).u{j} = vData.u{j}(index,:);
            data.(['p', num2str(i)]).v{j} = vData.v{j}(index,:);
            data.(['p', num2str(i)]).w{j} = vData.w{j}(index,:);
            data.(['p', num2str(i)]).uPrime{j} = vData.uPrime{j}(index,:);
            data.(['p', num2str(i)]).vPrime{j} = vData.vPrime{j}(index,:);
            data.(['p', num2str(i)]).wPrime{j} = vData.wPrime{j}(index,:);
        end
        
        % Revert Data Origin
        if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
            data.(['p', num2str(i)]).position(:,1) = data.(['p', num2str(i)]).position(:,1) - 1.325;

            if strcmp(data.(['p', num2str(i)]).planeOrientation, 'X')
                data.(['p', num2str(i)]).planePosition(:,1) = data.(['p', num2str(i)]).planePosition(:,1) -1.325;
            end
        end
        
        % Rename Plane
        pName = [data.(['p', num2str(i)]).planeOrientation, erase(num2str(data.(['p', num2str(i)]).planePosition, '%.5f'), '.')];
        data = renameStructField(data, ['p', num2str(i)], pName);
        data.(pName) = orderfields(data.(pName));
        
    end
    
    data = orderfields(data);
    
end


%% Local Functions

function pos = inputPos(type, plane)

    valid = false;
    while ~valid
        pos = str2double(input(['    ', type, ' ', plane, '-Position [m]: '], 's'));
        
        if isnan(pos) || length(pos) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end