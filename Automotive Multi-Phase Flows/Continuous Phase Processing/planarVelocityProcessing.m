%% Planar Velocity Processing v6.0
% ----
% Load, Process and Present Experimental and Numerical Planar Velocity Profiles


%% Preamble

run preamble;

%#ok<*UNRCH>

normDims = true; % Normalise Spatial Dimensions

figSave = false; % Save .fig File(s)

disp('===============================');
disp('Planar Velocity Processing v6.0');
disp('===============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Improve Versatility
% v3.0 - Moved Plotting to velocityPlots.m
% v4.0 - Rewrite to Support ParaView, Probe and Experimental Planar Data
% v5.0 - Update To Improve Consistency of Structures Across Repository
% v5.1 - Minor Update to Shift Preamble Into Separate Script
% v5.2 - Updates To Correct Inconsistent Normalisation Throughout Repository
% v5.3 - Removed Blockage Correction
% v6.0 - Added Support for Instantaneous Velocity Processing


%% Select Data Format

disp('Data Format');
disp('------------');

disp(' ');

disp('Compatible Data Formats:');
disp('    A: ParaView Planar Data');
disp('    B: Probe Planar Data');
disp('    C: Experimental Planar Data (Varney)');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Data Format [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';

        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';

        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Acquire Velocity Data

switch format

    case 'A'
        [campaignID, caseID, uData] = initialisePVdata(saveLoc, 'U');
        
        planes = fieldnames(uData);
        
    case 'B'
        maxNumCompThreads(nProc);
        
        [caseFolder, campaignID, caseID, timeDirs, deltaT, ...
         timePrecision, dataID, uData, sampleInt] = initialiseProbeData(saveLoc, maxNumCompThreads, ...
                                                                        'uProbes', 'wake');
        
        nTimes = height(uData.time);
        
    case 'C'
        [campaignID, caseID, uData] = initialiseExpFlowData(saveLoc, 'U');
        
        planes = fieldnames(uData);
        
        for i = 1:height(planes)
            uData.(planes{i}).u = rmfield(uData.(planes{i}).u, 'RMS');
            uData.(planes{i}).v = rmfield(uData.(planes{i}).v, 'RMS');
            uData.(planes{i}).w = rmfield(uData.(planes{i}).w, 'RMS');
        end
        
end

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');


%% Select Plane(s) of Interest

switch format
    
    case {'A', 'C'}
        
        valid = false;
        while ~valid
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'multiple', ...
                                     'name', 'Select Plane(s) to Process', ...
                                     'listString', planes);
        
            if ~valid
                disp(    'WARNING: No Planes Selected');
            end
        end
        clear valid;
        
        planes = planes(index);
        
    case 'B'
        uData = planarProbeDataExtraction(uData, spacePrecision, 'singlePlane');
        
        planes = fieldnames(uData);
        
end

disp(' ');
disp(' ');


%% Generate Velocity Profiles

disp('Velocity Profile Generation');
disp('----------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

%%%%

disp(' ');

% Adjust Data Origin
disp('    Adjusting Data Origin...');

switch format

    case {'A', 'B'}

        if strcmp(campaignID, 'Windsor_Upstream_2023')
            
            for i = 1:height(planes)
                uData.(planes{i}).positionGrid(:,1) = uData.(planes{i}).positionGrid(:,1) + 1.325;
            end
            clear i;

        end

end

disp(' ');

% Remove Duplicate Entries
disp('    Removing Duplicate Entries...');

for i = 1:height(planes)
    [uData.(planes{i}).positionGrid, index] = unique(uData.(planes{i}).positionGrid, 'rows', 'stable');
    
    switch format
        
        case {'A', 'C'}
            uData.(planes{i}).u.mean = uData.(planes{i}).u.mean(index);
            uData.(planes{i}).v.mean = uData.(planes{i}).v.mean(index);
            uData.(planes{i}).w.mean = uData.(planes{i}).w.mean(index);
        
        case 'B'
            
            for j = 1:nTimes
                uData.(planes{i}).u.inst{j} = uData.(planes{i}).u.inst{j}(index);
                uData.(planes{i}).v.inst{j} = uData.(planes{i}).v.inst{j}(index);
                uData.(planes{i}).w.inst{j} = uData.(planes{i}).w.inst{j}(index);
            end
            clear j;
            
    end
    
end
clear i;

disp(' ');

% Map Raw Data Onto Uniform Grid
disp('    Mapping Raw Data Onto Uniform Grid...');

% Set Target Spatial Resolution
if strcmp(campaignID, 'Windsor_fullScale')
    targetSize = 4e-3;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    targetSize = 1e-3;
else
    targetSize = 1e-3;
end

for i = 1:height(planes)

    if contains(planes{i}, 'X_')
        orientation = 'YZ';
    elseif contains(planes{i}, 'Y_')
        orientation = 'XZ';
    else
        orientation = 'XY';
    end

    switch orientation

        case 'YZ'
            xLimsData = double(uData.(planes{i}).positionGrid(1,1));
            yLimsData = double([min(uData.(planes{i}).positionGrid(:,2)); ...
                                max(uData.(planes{i}).positionGrid(:,2))]);
            zLimsData = double([min(uData.(planes{i}).positionGrid(:,3)); ...
                                max(uData.(planes{i}).positionGrid(:,3))]);
            
            % Adjust Uniform Cell Size to Fit Region of Interest
            nPy = (diff(yLimsData) / targetSize) + 1;
            nPz = (diff(zLimsData) / targetSize) + 1;

            sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
            sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));
            
            cellSize.(planes{i}).target = targetSize;
            cellSize.(planes{i}).x = targetSize;
            cellSize.(planes{i}).y = sizeY(1); clear sizeY;
            cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
            cellSize.(planes{i}).area = cellSize.(planes{i}).y * cellSize.(planes{i}).z;

            yOrig = double(uData.(planes{i}).positionGrid(:,2));
            zOrig = double(uData.(planes{i}).positionGrid(:,3));

            % Generate Grid
            [y, z] = ndgrid(linspace(yLimsData(1), yLimsData(2), nPy), ...
                            linspace(zLimsData(1), zLimsData(2), nPz));

            uData.(planes{i}).positionGrid = zeros([height(y(:)),3]);
            uData.(planes{i}).positionGrid(:,1) = xLimsData;
            uData.(planes{i}).positionGrid(:,[2,3]) = [y(:), z(:)];

            % Perform Interpolation
            switch format
                
                case {'A', 'C'}
                    uInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).u.mean, ...
                                                   'linear', 'none');
                    vInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).v.mean, ...
                                                   'linear', 'none');
                    wInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).w.mean, ...
                                                   'linear', 'none');

                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));

                case 'B'
                    
                    for j = 1:nTimes
                        uInterp = scatteredInterpolant(yOrig, zOrig, double(uData.(planes{i}).u.inst{j}), ...
                                                       'linear', 'none');
                        vInterp = scatteredInterpolant(yOrig, zOrig, double(uData.(planes{i}).v.inst{j}), ...
                                                       'linear', 'none');
                        wInterp = scatteredInterpolant(yOrig, zOrig, double(uData.(planes{i}).w.inst{j}), ...
                                                       'linear', 'none');

                        uData.(planes{i}).u.inst{j} = single(uInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                        uData.(planes{i}).v.inst{j} = single(vInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                        uData.(planes{i}).w.inst{j} = single(wInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                    end
                    clear j;
                    
            end

        case 'XZ'
            
            switch format
                
                case {'A', 'B'}
                    xLimsData = double([min(uData.(planes{i}).positionGrid(:,1)); ...
                                        max(uData.(planes{i}).positionGrid(:,1))]);
                    yLimsData = double(uData.(planes{i}).positionGrid(1,2));
                    zLimsData = double([min(uData.(planes{i}).positionGrid(:,3)); ...
                                        max(uData.(planes{i}).positionGrid(:,3))]);
                    
                case 'C'
                    xLimsData = [0.48525; 1.15];
                    yLimsData = uData.(planes{i}).positionGrid(1,2);
                    zLimsData = [0.035; 0.354];
                    
            end

            % Adjust Uniform Cell Size to Fit Region of Interest
            nPx = (diff(xLimsData) / targetSize) + 1;
            nPz = (diff(zLimsData) / targetSize) + 1;

            sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
            sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

            cellSize.(planes{i}).x = sizeX(1); clear sizeX;
            cellSize.(planes{i}).y = targetSize;
            cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
            cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).z;

            xOrig = double(uData.(planes{i}).positionGrid(:,1));
            zOrig = double(uData.(planes{i}).positionGrid(:,3));

            % Generate Grid
            [x, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                            linspace(zLimsData(1), zLimsData(2), nPz));

            uData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
            uData.(planes{i}).positionGrid(:,2) = yLimsData;
            uData.(planes{i}).positionGrid(:,[1,3]) = [x(:), z(:)];
            
            % Perform Interpolation
            switch format
                
                case {'A', 'C'}
                    uInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).u.mean, ...
                                                   'linear', 'none');
                    vInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).v.mean, ...
                                                   'linear', 'none');
                    wInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).w.mean, ...
                                                   'linear', 'none');

                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));

                case 'B'
                    
                    for j = 1:nTimes
                        uInterp = scatteredInterpolant(xOrig, zOrig, double(uData.(planes{i}).u.inst{j}), ...
                                                       'linear', 'none');
                        vInterp = scatteredInterpolant(xOrig, zOrig, double(uData.(planes{i}).v.inst{j}), ...
                                                       'linear', 'none');
                        wInterp = scatteredInterpolant(xOrig, zOrig, double(uData.(planes{i}).w.inst{j}), ...
                                                       'linear', 'none');

                        uData.(planes{i}).u.inst{j} = single(uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                        uData.(planes{i}).v.inst{j} = single(vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                        uData.(planes{i}).w.inst{j} = single(wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,3)));
                    end
                    clear j;
                    
            end

        case 'XY'

            switch format

                case {'A', 'B'}
                    xLimsData = double([min(uData.(planes{i}).positionGrid(:,1)); ...
                                        max(uData.(planes{i}).positionGrid(:,1))]);
                    yLimsData = double([min(uData.(planes{i}).positionGrid(:,2)); ...
                                        max(uData.(planes{i}).positionGrid(:,2))]);
                    zLimsData = double(uData.(planes{i}).positionGrid(1,3));
                    
                case 'C'
                    xLimsData = [0.48525; 1.15];
                    yLimsData = [-0.21; 0.21];
                    zLimsData = uData.(planes{i}).positionGrid(1,3);
                    
            end

            % Adjust Uniform Cell Size to Fit Region of Interest
            nPx = (diff(xLimsData) / targetSize) + 1;
            nPy = (diff(yLimsData) / targetSize) + 1;

            sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
            sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));

            cellSize.(planes{i}).x = sizeX(1); clear sizeX;
            cellSize.(planes{i}).y = sizeY(1); clear sizeY;
            cellSize.(planes{i}).z = targetSize;
            cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).y;

            xOrig = double(uData.(planes{i}).positionGrid(:,1));
            yOrig = double(uData.(planes{i}).positionGrid(:,2));

            % Generate Grid
            [x, y] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                            linspace(yLimsData(1), yLimsData(2), nPy));

            uData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
            uData.(planes{i}).positionGrid(:,3) = zLimsData;
            uData.(planes{i}).positionGrid(:,[1,2]) = [x(:), y(:)];
            
            % Perform Interpolation
            switch format
                
                case {'A', 'C'}
                    uInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).u.mean, ...
                                                   'linear', 'none');
                    vInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).v.mean, ...
                                                   'linear', 'none');
                    wInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).w.mean, ...
                                                   'linear', 'none');

                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));

                case 'B'
                    
                    for j = 1:nTimes
                        uInterp = scatteredInterpolant(xOrig, yOrig, double(uData.(planes{i}).u.inst{j}), ...
                                                       'linear', 'none');
                        vInterp = scatteredInterpolant(xOrig, yOrig, double(uData.(planes{i}).v.inst{j}), ...
                                                       'linear', 'none');
                        wInterp = scatteredInterpolant(xOrig, yOrig, double(uData.(planes{i}).w.inst{j}), ...
                                                       'linear', 'none');

                        uData.(planes{i}).u.inst{j} = single(uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,2)));
                        uData.(planes{i}).v.inst{j} = single(vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,2)));
                        uData.(planes{i}).w.inst{j} = single(wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                                     uData.(planes{i}).positionGrid(:,2)));
                    end
                    clear j;
                    
            end

    end
    clear targetSize orientation xProg yOrig zOrig x y z uInterp vInterp wInterp;

end
clear i;

% Calculate Time-Averaged and Fluctuating Fields
switch format
    
    case 'B'
        disp(' ');
        
        disp('    Calculating Time-Averaged Field(s)...');
        
        for i = 1:height(planes)
            uData.(planes{i}).u.mean = zeros([height(uData.(planes{i}).positionGrid),1], 'single');
            uData.(planes{i}).v.mean = uData.(planes{i}).u.mean;
            uData.(planes{i}).w.mean = uData.(planes{i}).u.mean;
            
            for j = 1:nTimes
                uData.(planes{i}).u.mean = uData.(planes{i}).u.mean + uData.(planes{i}).u.inst{j};
                uData.(planes{i}).v.mean = uData.(planes{i}).v.mean + uData.(planes{i}).v.inst{j};
                uData.(planes{i}).w.mean = uData.(planes{i}).w.mean + uData.(planes{i}).w.inst{j};
            end
            clear j;
            
            uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / nTimes;
            uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / nTimes;
            uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / nTimes;
        end
        clear i;
        
        disp(' ');
        
        disp('    Calculating Instantaneous Field Fluctuations...');
        
        for i = 1:height(planes)
            uData.(planes{i}).u.prime = uData.(planes{i}).u.inst;
            uData.(planes{i}).v.prime = uData.(planes{i}).v.inst;
            uData.(planes{i}).w.prime = uData.(planes{i}).w.inst;
            
            for j = 1:nTimes
                uData.(planes{i}).u.prime{j} = uData.(planes{i}).u.prime{j} - uData.(planes{i}).u.mean;
                uData.(planes{i}).v.prime{j} = uData.(planes{i}).v.prime{j} - uData.(planes{i}).v.mean;
                uData.(planes{i}).w.prime{j} = uData.(planes{i}).w.prime{j} - uData.(planes{i}).w.mean;
            end
            clear j;
            
        end
        clear i;
        
        disp(' ');
        
        disp('    Calculating RMS of Field Fluctuations...');
        
        for i = 1:height(planes)
            uData.(planes{i}).u.RMS = zeros([height(uData.(planes{i}).positionGrid),1], 'single');
            uData.(planes{i}).v.RMS = uData.(planes{i}).u.RMS;
            uData.(planes{i}).w.RMS = uData.(planes{i}).u.RMS;
            
            for j = 1:nTimes
                uData.(planes{i}).u.RMS = uData.(planes{i}).u.RMS + uData.(planes{i}).u.prime{j}.^2;
                uData.(planes{i}).v.RMS = uData.(planes{i}).v.RMS + uData.(planes{i}).v.prime{j}.^2;
                uData.(planes{i}).w.RMS = uData.(planes{i}).w.RMS + uData.(planes{i}).w.prime{j}.^2;
            end
            clear j;
            
            uData.(planes{i}).u.RMS = sqrt((1 / nTimes) * uData.(planes{i}).u.RMS);
            uData.(planes{i}).v.RMS = sqrt((1 / nTimes) * uData.(planes{i}).v.RMS);
            uData.(planes{i}).w.RMS = sqrt((1 / nTimes) * uData.(planes{i}).w.RMS);
        end
        clear i;

end

disp(' ');

% Calculate Vorticity
disp('    Calculating Vorticity...');

for i = 1:height(planes)

    if contains(planes{i}, 'X_')
        orientation = 'YZ';
    elseif contains(planes{i}, 'Y_')
        orientation = 'XZ';
    else
        orientation = 'XY';
    end

    switch orientation

        case 'YZ'
            gridShape = [height(unique(uData.(planes{i}).positionGrid(:,2))), ...
                         height(unique(uData.(planes{i}).positionGrid(:,3)))];

            y = reshape(uData.(planes{i}).positionGrid(:,2), gridShape)';
            z = reshape(uData.(planes{i}).positionGrid(:,3), gridShape)';

            v = reshape(uData.(planes{i}).v.mean, gridShape)';
            w = reshape(uData.(planes{i}).w.mean, gridShape)';

            [uData.(planes{i}).omega.mean, ~] = curl(y, z, v, w);
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean';
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean(:);

            % Omit Values Close to the Ground Plane
            if strcmp(campaignID, 'Windsor_fullScale')
                uData.(planes{i}).omega.mean(uData.(planes{i}).positionGrid(:,3) < 32e-3) = 0;
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                uData.(planes{i}).omega.mean(uData.(planes{i}).positionGrid(:,3) < 8e-3) = 0;
            end

        case 'XZ'
            gridShape = [height(unique(uData.(planes{i}).positionGrid(:,1))), ...
                         height(unique(uData.(planes{i}).positionGrid(:,3)))];

            x = reshape(uData.(planes{i}).positionGrid(:,1), gridShape)';
            z = reshape(uData.(planes{i}).positionGrid(:,3), gridShape)';

            u = reshape(uData.(planes{i}).u.mean, gridShape)';
            w = reshape(uData.(planes{i}).w.mean, gridShape)';

            [uData.(planes{i}).omega.mean, ~] = curl(x, z, u, w);
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean';
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean(:);

            % Omit Values Close to the Ground Plane
            if strcmp(campaignID, 'Windsor_fullScale')
                uData.(planes{i}).omega.mean(uData.(planes{i}).positionGrid(:,3) < 32e-3) = 0;
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                uData.(planes{i}).omega.mean(uData.(planes{i}).positionGrid(:,3) < 8e-3) = 0;
            end

        case 'XY'
            gridShape = [height(unique(uData.(planes{i}).positionGrid(:,1))), ...
                         height(unique(uData.(planes{i}).positionGrid(:,2)))];

            x = reshape(uData.(planes{i}).positionGrid(:,1), gridShape)';
            y = reshape(uData.(planes{i}).positionGrid(:,2), gridShape)';

            u = reshape(uData.(planes{i}).u.mean, gridShape)';
            v = reshape(uData.(planes{i}).v.mean, gridShape)';

            [uData.(planes{i}).omega.mean, ~] = curl(x, y, u, v);
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean';
            uData.(planes{i}).omega.mean = uData.(planes{i}).omega.mean(:);

    end

    % Omit Artificially High Values
    vortLims = prctile(uData.(planes{i}).omega.mean, [1, 99]);
    uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean < vortLims(1)) = vortLims(1);
    uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean > vortLims(2)) = vortLims(2);

    % Normalise
    uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean < 0) = ...
    rescale(uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean < 0), -1, 0);
    uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean > 0) = ...
    rescale(uData.(planes{i}).omega.mean(uData.(planes{i}).omega.mean > 0), 0, 1);
    uData.(planes{i}).omega.mean((uData.(planes{i}).omega.mean > -0.25) & ...
    (uData.(planes{i}).omega.mean < 0.25)) = 0;

    clear orientation x y z u v w vortLims;

end
clear i;

%%%%

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Save Instantaneous Data

switch format
    
    case 'B'
        disp('Data Save Options');
        disp('------------------');

        valid = false;
        while ~valid
            disp(' ');

            selection = input('Save Data for Future Use? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>

                % Save Data
                if ~exist([saveLoc, '/Numerical/MATLAB/planarVelocityFields/', campaignID, '/', caseID, '/', dataID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarVelocityFields/', campaignID, '/', caseID, '/', dataID]);
                end
                
                planeID = planes{1};
                
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarVelocityFields/', campaignID, '/', caseID, '/', dataID,  '/', planeID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarVelocityFields/', campaignID, '/', caseID, '/', dataID, '/', planeID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'planeID', 'uData', 'cellSize', 'sampleInt', 'timePrecision', '-v7.3', '-noCompression');
                disp('        Success');

                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;

        disp(' ');
        disp(' ');
        
end


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Field(s)? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMean = false;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMean = true;

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

switch format
    
    case 'B'
        valid = false;
        while ~valid
            disp(' ');

            selection = input('Plot RMS Field(s)? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotRMS = false;

                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                plotRMS = true;

                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;

        valid = false;
        while ~valid
            disp(' ');

            selection = input('Plot Instantaneous Fields? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotInst = false;

                valid = true;
            elseif selection == 'y' | selection == 'Y' %#ok<OR2>
                plotInst = true;

                startFrame = inputFrames(nTimes, 'Start');

                if startFrame == -1
                    continue;
                end

                endFrame = inputFrames(nTimes, 'End');

                if endFrame == -1
                    continue;
                elseif endFrame < startFrame
                    disp('        WARNING: Invalid Time Format (''endFrame'' Precedes ''startFrame'')');

                    continue;
                end

                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;
        
    otherwise
        plotRMS = false;
        plotInst = false;

end

if plotMean || plotInst
    
    valid = false;
    while ~valid
        disp(' ');

        selection = input('Plot Vorticity Field(s)? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            plotOmega = false;

            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            plotOmega = true;

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;
    
else
    plotRMS = false;
end

if plotMean || plotRMS || plotInst
    disp(' ');
    
    disp('Preparing Plots...')
    
    % Normalise Coordinate System
    if normDims
        disp('    Normalising Spatial Dimensions');

        parts = fieldnames(geometry);
        for i = 1:height(parts)
            geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
        end
        clear i parts;

        xDims = xDims / normLength;
        yDims = yDims / normLength;
        zDims = zDims / normLength;

        for i = 1:height(planes)
            cellSize.(planes{i}).target = cellSize.(planes{i}).target / normLength;
            cellSize.(planes{i}).x = cellSize.(planes{i}).x / normLength;
            cellSize.(planes{i}).y = cellSize.(planes{i}).y / normLength;
            cellSize.(planes{i}).z = cellSize.(planes{i}).z / normLength;
            cellSize.(planes{i}).area = cellSize.(planes{i}).area / (normLength^2);

            uData.(planes{i}).positionGrid = uData.(planes{i}).positionGrid / normLength;
        end

    end
    
    % Normalise Velocity
    disp('    Normalising Velocity');
    
    switch format
        
        case 'A'
            
            if strcmp(campaignID, 'Windsor_fullScale')
                U = 22.22; % m/s
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                U = 40; % m/s
            else
                U = 1;
            end
            
            for i = 1:height(planes)
                uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
                uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
                uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
            end
            clear i;
            
        case 'B'

            if strcmp(campaignID, 'Windsor_fullScale')
                U = 22.22; % m/s
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                U = 40; % m/s
            else
                U = 1;
            end
            
            for i = 1:height(planes)
                uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
                uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
                uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
                
                uData.(planes{i}).u.RMS = uData.(planes{i}).u.RMS / U;
                uData.(planes{i}).v.RMS = uData.(planes{i}).v.RMS / U;
                uData.(planes{i}).w.RMS = uData.(planes{i}).w.RMS / U;
                
                for j = 1:nTimes
                    uData.(planes{i}).u.inst{j} = uData.(planes{i}).u.inst{j} / U;
                    uData.(planes{i}).v.inst{j} = uData.(planes{i}).v.inst{j} / U;
                    uData.(planes{i}).w.inst{j} = uData.(planes{i}).w.inst{j} / U;
                    
                    uData.(planes{i}).u.prime{j} = uData.(planes{i}).u.prime{j} / U;
                    uData.(planes{i}).v.prime{j} = uData.(planes{i}).v.prime{j} / U;
                    uData.(planes{i}).w.prime{j} = uData.(planes{i}).w.prime{j} / U;
                end
                clear j;
                
            end
            clear i;

        case 'C'

            if strcmp(campaignID, 'Varney')
                U = 40; % m/s

                for i = 1:height(planes)
                        uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
                        uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
                        uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
                end

            end

    end
    
    if strcmp(campaignID, 'Windsor_fullScale')
        spatialRes = 2e-3;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        spatialRes = 0.5e-3;
    else
        spatialRes = 0.5e-3;
    end

    if normDims
        spatialRes = spatialRes / normLength;
    end

    component = [];
    mapPerim = [];
    nPlanes = 1;
    planeNo = 1;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
end

disp(' ');
disp(' ');


%% Present Velocity Data

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean

    for i = 1:height(planes)
        
        switch format
            
            case {'A', 'C'}
                disp(['Presenting ', planes{i}, '...']);
                
            case 'B'
                disp('Presenting the Time-Averaged Flow Field');
                
        end

        if contains(planes{i}, 'X_')
            orientation = 'YZ';
        elseif contains(planes{i}, 'Y_')
            orientation = 'XZ';
        else
            orientation = 'XY';
        end

        positionData = uData.(planes{i}).positionGrid;
        vectorData = [uData.(planes{i}).u.mean, uData.(planes{i}).v.mean, uData.(planes{i}).w.mean];

        switch orientation
            
            case 'YZ'
                xLimsPlot = [0.3; 4.625766283524905];
                yLimsPlot = [-0.5; 0.5];
                zLimsPlot = [0; 0.5];
                
            case {'XZ', 'XY'}
                xLimsPlot = [0.3; 1.2];
                yLimsPlot = [-0.3; 0.3];
                zLimsPlot = [0; 0.5];
                
%                 xLimsPlot = [0.3; 0.9];
%                 yLimsPlot = [0; 0.3];
%                 zLimsPlot = [0; 0.5];
                
        end

        if ~normDims
            xLimsPlot = xLimsPlot * normLength;
            yLimsPlot = yLimsPlot * normLength;
            zLimsPlot = zLimsPlot * normLength;
        end

        switch orientation

            case 'YZ'
                xLimsData = uData.(planes{i}).positionGrid(1,1);
                yLimsData = yLimsPlot;
                zLimsData = zLimsPlot;

            case 'XZ'
                xLimsData = xLimsPlot;
                yLimsData = uData.(planes{i}).positionGrid(1,2);
                zLimsData = zLimsPlot;

            case 'XY'
                xLimsData = xLimsPlot;
                yLimsData = yLimsPlot;
                zLimsData = uData.(planes{i}).positionGrid(1,3);

        end

        switch orientation

            case 'YZ'
                nComponents = 3;

            case {'XZ', 'XY'}
                nComponents = 2;

        end
        
        switch format
            
            case {'A', 'C'}
                figName = [planes{i}, '_', caseID];
                
            case 'B'
                figName = [planes{i}, '_U_Mean_', caseID];
                
        end
                
        cMap = viridis(32);
        streamlines = true;
        cLims = [0, 1];

        [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                               streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);

        if plotOmega
            scalarData = uData.(planes{i}).omega.mean;
            figName = [planes{i}, '_Vorticity_', caseID];
            cMap = cool2warm(32);
            contourlines = [];
            refPoint = [];
            cLims = [-1; 1];

            [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                   planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                   refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                   zLimsPlot, normDims, figSave);
        end
        
    end
    clear i;
    
    disp(' ');    
end

if plotRMS
    disp('    Presenting the RMS of the Fluctuating Flow Field');
    
    for i = 1:height(planes)
        
        if contains(planes{i}, 'X_')
            orientation = 'YZ';
        elseif contains(planes{i}, 'Y_')
            orientation = 'XZ';
        else
            orientation = 'XY';
        end

        positionData = uData.(planes{i}).positionGrid;
        scalarData = sqrt(uData.(planes{i}).u.RMS.^2 + ...
                          uData.(planes{i}).v.RMS.^2 + ...
                          uData.(planes{i}).w.RMS.^2);

        switch orientation
            
            case 'YZ'
                xLimsPlot = [0.3; 4.625766283524905];
                yLimsPlot = [-0.5; 0.5];
                zLimsPlot = [0; 0.5];
                
            case {'XZ', 'XY'}
                xLimsPlot = [0.3; 1.2];
                yLimsPlot = [-0.3; 0.3];
                zLimsPlot = [0; 0.5];
                
%                 xLimsPlot = [0.3; 0.9];
%                 yLimsPlot = [0; 0.3];
%                 zLimsPlot = [0; 0.5];
                
        end

        if ~normDims
            xLimsPlot = xLimsPlot * normLength;
            yLimsPlot = yLimsPlot * normLength;
            zLimsPlot = zLimsPlot * normLength;
        end

        switch orientation

            case 'YZ'
                xLimsData = uData.(planes{i}).positionGrid(1,1);
                yLimsData = yLimsPlot;
                zLimsData = zLimsPlot;

            case 'XZ'
                xLimsData = xLimsPlot;
                yLimsData = uData.(planes{i}).positionGrid(1,2);
                zLimsData = zLimsPlot;

            case 'XY'
                xLimsData = xLimsPlot;
                yLimsData = yLimsPlot;
                zLimsData = uData.(planes{i}).positionGrid(1,3);

        end
        
        figName = [planes{i}, '_RMS_U_', caseID];
        cMap = viridis(32);
        contourlines = [];
        refPoint = [];
        cLims = [0; max(scalarData)];

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear i;
    
    disp(' ');    
end

if plotInst

    for i = 1:height(planes)
        disp('    Presenting the Instantaneous Flow Field');

        if contains(planes{i}, 'X_')
            orientation = 'YZ';
        elseif contains(planes{i}, 'Y_')
            orientation = 'XZ';
        else
            orientation = 'XY';
        end

        positionData = uData.(planes{i}).positionGrid;

        switch orientation
            
            case 'YZ'
                xLimsPlot = [0.3; 4.625766283524905];
                yLimsPlot = [-0.5; 0.5];
                zLimsPlot = [0; 0.5];
                
            case {'XZ', 'XY'}
                xLimsPlot = [0.3; 1.2];
                yLimsPlot = [-0.3; 0.3];
                zLimsPlot = [0; 0.5];
                
%                 xLimsPlot = [0.3; 0.9];
%                 yLimsPlot = [0; 0.3];
%                 zLimsPlot = [0; 0.5];
                
        end

        if ~normDims
            xLimsPlot = xLimsPlot * normLength;
            yLimsPlot = yLimsPlot * normLength;
            zLimsPlot = zLimsPlot * normLength;
        end

        switch orientation

            case 'YZ'
                xLimsData = uData.(planes{i}).positionGrid(1,1);
                yLimsData = yLimsPlot;
                zLimsData = zLimsPlot;

            case 'XZ'
                xLimsData = xLimsPlot;
                yLimsData = uData.(planes{i}).positionGrid(1,2);
                zLimsData = zLimsPlot;

            case 'XY'
                xLimsData = xLimsPlot;
                yLimsData = yLimsPlot;
                zLimsData = uData.(planes{i}).positionGrid(1,3);

        end

        switch orientation

            case 'YZ'
                nComponents = 3;

            case {'XZ', 'XY'}
                nComponents = 2;

        end
                
        cMap = viridis(32);
        streamlines = true;
        cLims = [0, 1];
        
        figHold = fig;
        
        for j = startFrame:endFrame
            
            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            vectorData = [uData.(planes{i}).u.inst{j}, ...
                          uData.(planes{i}).v.inst{j}, ...
                          uData.(planes{i}).w.inst{j}];
            figTime = num2str(uData.(planes{i}).time(j), ['%.', num2str(timePrecision), 'f']);
            figName = [planes{i}, '_Inst_U_', '_T', erase(figTime, '.'), '_', caseID];
            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                                   mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                   streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                   zLimsPlot, normDims, figSave);  
        end
        clear j;
        
    end
    clear i;
    
end

if ~plotMean && ~ plotRMS && ~plotInst
    disp('Skipping Map Presentation...');
end

disp(' ');
disp(' ');


%% Local Functions

function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
