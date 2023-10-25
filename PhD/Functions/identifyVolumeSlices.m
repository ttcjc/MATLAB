%% Volume Slice Identifier v1.0
% ----
% Extract 2D Planes From an Existing 3D Grid
% ----
% Usage: volumeSlice = identifyVolumeSlices(positionData, spacePrecision, multiSlice)
% 
%        'positionData'   -> 
%        'spacePrecision' -> 
%        'multiSlice'     -> 


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function volumeSlice = identifyVolumeSlices(positionData, spacePrecision, multiSlice)

    disp('Plane of Interest Selection');
    disp('----------------------------');
    
    disp(' ');
    
    % Identify Data Limits
    xLimsData = [min(positionData(:,1)); max(positionData(:,1))];
    yLimsData = [min(positionData(:,2)); max(positionData(:,2))];
    zLimsData = [min(positionData(:,3)); max(positionData(:,3))];
    
    disp('Volume Boundaries:');
    disp(['    X: ', num2str(xLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(xLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
    disp(['    Y: ', num2str(yLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(yLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
    disp(['    Z: ', num2str(zLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(zLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);

    disp(' ');

    disp('Possible Planar Orientations:');
    disp('    X (YZ-Plane)');
    disp('    Y (XZ-Plane)');
    disp('    Z (XY-Plane)');
    
    % Specify Slice Count
    if multiSlice
        
        valid = false;
        while ~valid
            disp(' ');
            nSlices = str2num(input('Select Desired Number of Slices [Scalar]: ', 's')); %#ok<ST2NM>
            
            if isempty(nSlices) || isnan(nSlices) || length(nSlices) > 1
                disp('    WARNING: Invalid Entry');
            else
                valid = true;
            end
        end
        
    else
        nSlices = 1;
    end
    
    for i = 1:nSlices
        
        % Define Slice Orientation
        valid = false;
        while ~valid
            disp(' ');
            
            if nSlices == 1
                selection = input('Select Orientation of Slice [X/Y/Z]: ', 's');
            else
                selection = input(['Select Orientation of Slice #', num2str(i), ' [X/Y/Z]: '], 's');
            end

            if selection == 'x' | selection == 'X' %#ok<OR2>
                volumeSlice.(['s', num2str(i)]).orientation = 'YZ';
                
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                volumeSlice.(['s', num2str(i)]).orientation = 'XZ';
                
                valid = true;
            elseif selection == 'z' | selection == 'Z' %#ok<OR2>
                volumeSlice.(['s', num2str(i)]).orientation = 'XY';
                
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
        volumeSlice.(['s', num2str(i)]).xLims = zeros(2,1);
        volumeSlice.(['s', num2str(i)]).yLims = volumeSlice.(['s', num2str(i)]).xLims;
        volumeSlice.(['s', num2str(i)]).zLims = volumeSlice.(['s', num2str(i)]).xLims;

        % Define Slice Limits
        valid = false;
        while ~valid
            
            switch volumeSlice.(['s', num2str(i)]).orientation
                
                case 'YZ'                    
                    volumeSlice.(['s', num2str(i)]).xLims = inputPos('Planar', 'X');
                    
                    if volumeSlice.(['s', num2str(i)]).xLims == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    
                    if volumeSlice.(['s', num2str(i)]).yLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    
                    if volumeSlice.(['s', num2str(i)]).yLims(2) == -1
                        continue;
                    end         
                    
                    volumeSlice.(['s', num2str(i)]).yLims = sort(volumeSlice.(['s', num2str(i)]).yLims);
                    
                    volumeSlice.(['s', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    
                    if volumeSlice.(['s', num2str(i)]).zLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    
                    if volumeSlice.(['s', num2str(i)]).zLims(2) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).zLims = sort(volumeSlice.(['s', num2str(i)]).zLims);
                
                case 'XZ'                    
                    volumeSlice.(['s', num2str(i)]).yLims = inputPos('Planar', 'Y');
                    
                    if volumeSlice.(['s', num2str(i)]).yLims == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    
                    if volumeSlice.(['s', num2str(i)]).xLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    
                    if volumeSlice.(['s', num2str(i)]).xLims(2) == -1
                        continue;
                    end         
                    
                    volumeSlice.(['s', num2str(i)]).xLims = sort(volumeSlice.(['s', num2str(i)]).xLims);
                    
                    volumeSlice.(['s', num2str(i)]).zLims(1) = inputPos('Lower', 'Z');
                    
                    if volumeSlice.(['s', num2str(i)]).zLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).zLims(2) = inputPos('Upper', 'Z');
                    
                    if volumeSlice.(['s', num2str(i)]).zLims(2) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).zLims = sort(volumeSlice.(['s', num2str(i)]).zLims);
                    
                case 'XY'
                    volumeSlice.(['s', num2str(i)]).zLims = inputPos('Planar', 'Z');
                    
                    if volumeSlice.(['s', num2str(i)]).zLims == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).xLims(1) = inputPos('Lower', 'X');
                    
                    if volumeSlice.(['s', num2str(i)]).xLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).xLims(2) = inputPos('Upper', 'X');
                    
                    if volumeSlice.(['s', num2str(i)]).xLims(2) == -1
                        continue;
                    end         
                    
                    volumeSlice.(['s', num2str(i)]).xLims = sort(volumeSlice.(['s', num2str(i)]).xLims);
                    
                    volumeSlice.(['s', num2str(i)]).yLims(1) = inputPos('Lower', 'Y');
                    
                    if volumeSlice.(['s', num2str(i)]).yLims(1) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).yLims(2) = inputPos('Upper', 'Y');
                    
                    if volumeSlice.(['s', num2str(i)]).yLims(2) == -1
                        continue;
                    end
                    
                    volumeSlice.(['s', num2str(i)]).yLims = sort(volumeSlice.(['s', num2str(i)]).yLims);
            
            end
            
            if (min(volumeSlice.(['s', num2str(i)]).xLims) < xLimsData(1)) || (max(volumeSlice.(['s', num2str(i)]).xLims) > xLimsData(2)) || ...
               (min(volumeSlice.(['s', num2str(i)]).yLims) < yLimsData(1)) || (max(volumeSlice.(['s', num2str(i)]).yLims) > yLimsData(2)) || ...
               (min(volumeSlice.(['s', num2str(i)]).zLims) < zLimsData(1)) || (max(volumeSlice.(['s', num2str(i)]).zLims) > zLimsData(2))
                disp('        WARNING: Plane Lies Outside Volume');
                disp(' ');
            else
                valid = true;
            end
            
        end

        % Shift Requested Slice to Nearest Data Point
        switch volumeSlice.(['s', num2str(i)]).orientation
            
            case 'YZ'
                allPlanes = unique(positionData(:,1), 'stable');
                [offset, index] = min(abs(allPlanes - volumeSlice.(['s', num2str(i)]).xLims));
                
                if offset ~= 0
                    disp(' ');
                    disp('    INFO: Requested Slice Unavailable');
                    disp(['        Shifting X: ', num2str(volumeSlice.(['s', num2str(i)]).xLims, ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', ...
                                                  num2str(allPlanes(index), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
                    
                    volumeSlice.(['s', num2str(i)]).xLims = allPlanes(index);
                end
                
                volumeSlice.(['s', num2str(i)]).position = volumeSlice.(['s', num2str(i)]).xLims;
            
            case 'XZ'
                allPlanes = unique(positionData(:,2), 'stable');
                [offset, index] = min(abs(allPlanes - volumeSlice.(['s', num2str(i)]).yLims));
                
                if offset ~= 0
                    disp(' ');
                    disp('    INFO: Requested Slice Unavailable');
                    disp(['        Shifting Y: ', num2str(volumeSlice.(['s', num2str(i)]).yLims, ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', ...
                                                  num2str(allPlanes(index), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
                    
                    volumeSlice.(['s', num2str(i)]).yLims = allPlanes(index);
                end
                
                volumeSlice.(['s', num2str(i)]).position = volumeSlice.(['s', num2str(i)]).yLims;
            
            case 'XY'
                allPlanes = unique(positionData(:,3), 'stable');
                [offset, index] = min(abs(allPlanes - volumeSlice.(['s', num2str(i)]).zLims));
                
                if offset ~= 0
                    disp(' ');
                    disp('    INFO: Requested Slice Unavailable');
                    disp(['        Shifting Z: ', num2str(volumeSlice.(['s', num2str(i)]).zLims, ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', ...
                                                  num2str(allPlanes(index), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
                    
                    volumeSlice.(['s', num2str(i)]).zLims = allPlanes(index);
                end
                
                volumeSlice.(['s', num2str(i)]).position = volumeSlice.(['s', num2str(i)]).zLims;
        
        end

        % Rename Slice
        if volumeSlice.(['s', num2str(i)]).position > 0
            sliceName = [volumeSlice.(['s', num2str(i)]).orientation, '_P', erase(num2str(volumeSlice.(['s', num2str(i)]).position, ['%.', num2str(spacePrecision), 'f']), '.')];
        elseif volumeSlice.(['s', num2str(i)]).position < 0
            sliceName = [volumeSlice.(['s', num2str(i)]).orientation, '_N', erase(num2str(abs(volumeSlice.(['s', num2str(i)]).position), ['%.', num2str(spacePrecision), 'f']), '.')];
        else
            sliceName = [volumeSlice.(['s', num2str(i)]).orientation, '_0'];
        end

        volumeSlice = renameStructField(volumeSlice, ['s', num2str(i)], sliceName);
        volumeSlice.(sliceName) = orderfields(volumeSlice.(sliceName), {'orientation', 'position', 'xLims', 'yLims', 'zLims'});
        
    end
    
    if ~multiSlice
        volumeSlice = volumeSlice.(sliceName);
    end

end


%% Local Functions

function pos = inputPos(type, plane)

    pos = str2double(input(['    ', type, ' ', plane, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        
        pos = -1;
    end
    
end