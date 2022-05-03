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

function vData = extractPlanarProbeData(caseFolder, vData)

    % Temporarily Adjust Data Origin
    if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
        vData.position(:,1) = vData.position(:,1) + 1.325;
    end
    
    % Identify Data Limits
    xLims = [min(vData.position(:,1)); max(vData.position(:,1))];
    yLims = [min(vData.position(:,2)); max(vData.position(:,2))];
    zLims = [min(vData.position(:,3)); max(vData.position(:,3))];

    disp('Plane Specification');
    disp('--------------------');

    disp(' ');
    
    disp('Available Probe Volume:');
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
        nPlanes = str2num(input('Select Desired Number of Planes (Numerical Format): ', 's')); %#ok<ST2NM>

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
            selection = input('Select Orientation of Plane ', num2str(i), ' [X/Y/Z]: ', 's');

            if selection == 'x' | selection == 'X' %#ok<OR2>
                data.(i).planeOrientation = 'X';
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                data.(i).planeOrientation = 'Y';
                valid = true;
            elseif selection == 'z' | selection == 'Z' %#ok<OR2>
                data.(i).planeOrientation = 'Z';
                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        
    end
    
end


%% Local Functions