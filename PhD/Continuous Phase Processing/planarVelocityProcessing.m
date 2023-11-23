%% Planar Velocity Processing v5.2
% ----
% Load, Process and Present Experimental and Numerical Planar Velocity Profiles


%% Preamble

run preamble;

%#ok<*UNRCH>

normDims = true; % Normalise Spatial Dimensions

figSave = false; % Save .fig File(s)

disp('===============================');
disp('Planar Velocity Processing v5.2');
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
        error('NYI');
        
    case 'C'
        [campaignID, caseID, uData] = initialiseExpFlowData(saveLoc, 'U');
        
        planes = fieldnames(uData);
        
        for i = 1:height(planes)
            uData.(planes{i}).u = rmfield(uData.(planes{i}).u, 'RMS');
            uData.(planes{i}).v = rmfield(uData.(planes{i}).v, 'RMS');
            uData.(planes{i}).w = rmfield(uData.(planes{i}).w, 'RMS');
        end
        
end

%     case 'B'
%         [caseName, dataID, uData, sampleInterval, timePrecision, geometry, ...
%          xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(saveLoc, 'planar', normalise, nProc);    

%     case 'C'
%         [caseName, uData, geometry, ...
%          xDims, yDims, zDims, spacePrecision] = initialiseExpData('U', normalise);

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry;

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

    case 'A'

        if strcmp(campaignID, 'Windsor_Upstream_2023')
            
            for i = 1:height(planes)
                uData.(planes{i}).positionGrid(:,1) = uData.(planes{i}).positionGrid(:,1) + 1.325;
            end
            clear i;

        end

end

%     case 'B'
%         
%         if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
%             
%             for i = 1:height(planes)
%                 
%                 if normalise
%                     uData.(planes{i}).position(:,1) = uData.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);
%                     
%                     switch uData.(planes{i}).planeOrientation
%                         
%                         case 'YZ'
%                             uData.(planes{i}).planePosition = uData.(planes{i}).planePosition + round((1.325 / 1.044), spacePrecision);
%                             uData.(planes{i}).xLims = uData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
%                             
%                         case {'XZ', 'XY'}
%                             uData.(planes{i}).xLims = uData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
%                             
%                     end
%                     
%                 else
%                     uData.(planes{i}).position(:,1) = uData.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>
%                     
%                     switch uData.(planes{i}).planeOrientation
%                         
%                         case 'YZ'
%                             uData.(planes{i}).planePosition = uData.(planes{i}).planePosition + 1.325;
%                             uData.(planes{i}).xLims = uData.(planes{i}).xLims + 1.325;
%                         
%                         case {'XZ', 'XY'}
%                             uData.(planes{i}).xLims = uData.(planes{i}).xLims + 1.325;
%                     
%                     end
%                     
%                 end
% 
%             end
%             
%

for i = 1:height(planes)
    [uData.(planes{i}).positionGrid, index] = unique(uData.(planes{i}).positionGrid, 'rows', 'stable');

    uData.(planes{i}).u.mean = uData.(planes{i}).u.mean(index);
    uData.(planes{i}).v.mean = uData.(planes{i}).v.mean(index);
    uData.(planes{i}).w.mean = uData.(planes{i}).w.mean(index);
end
clear i;

disp(' ');

% Map Raw Data Onto UniformGrid
disp('    Mapping Raw Data Onto Uniform Grid...');

if strcmp(campaignID, 'Windsor_fullScale')
    cellSize.target = 4e-3;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    cellSize.target = 1e-3;
else
    cellSize.target = 1e-3;
end

switch format

    case {'A', 'C'}

        for i = 1:height(planes)
        
            if contains(planes{i}, '_X_')
                orientation = 'YZ';
            elseif contains(planes{i}, '_Y_')
                orientation = 'XZ';
            else
                orientation = 'XY';
            end
        
            switch orientation
        
                case 'YZ'
                    xLimsData = uData.(planes{i}).positionGrid(1,1);
                    yLimsData = [min(uData.(planes{i}).positionGrid(:,2)); ...
                                 max(uData.(planes{i}).positionGrid(:,2))];
                    zLimsData = [min(uData.(planes{i}).positionGrid(:,3)); ...
                                 max(uData.(planes{i}).positionGrid(:,3))];

                    nPy = (diff(yLimsData) / cellSize.target) + 1;
                    nPz = (diff(zLimsData) / cellSize.target) + 1;
                    
                    sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
                    sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

                    cellSize.(planes{i}).x = cellSize.target;
                    cellSize.(planes{i}).y = sizeY(1); clear sizeY;
                    cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
                    cellSize.(planes{i}).area = cellSize.(planes{i}).y * cellSize.(planes{i}).z;
        
                    yOrig = uData.(planes{i}).positionGrid(:,2);
                    zOrig = uData.(planes{i}).positionGrid(:,3);
                    
                    [y, z] = ndgrid(linspace(yLimsData(1), yLimsData(2), nPy), ...
                                    linspace(zLimsData(1), zLimsData(2), nPz));
                    
                    uData.(planes{i}).positionGrid = zeros([height(y(:)),3]);
                    uData.(planes{i}).positionGrid(:,1) = xLimsData;
                    uData.(planes{i}).positionGrid(:,[2,3]) = [y(:), z(:)];
                    
                    uInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).u.mean, 'linear', 'none');
                    vInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).v.mean, 'linear', 'none');
                    wInterp = scatteredInterpolant(yOrig, zOrig, uData.(planes{i}).w.mean, 'linear', 'none');
                    
                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,2), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    
                case 'XZ'
                    xLimsData = [min(uData.(planes{i}).positionGrid(:,1)); ...
                                 max(uData.(planes{i}).positionGrid(:,1))];
                    yLimsData = uData.(planes{i}).positionGrid(1,2);
                    zLimsData = [min(uData.(planes{i}).positionGrid(:,3)); ...
                                 max(uData.(planes{i}).positionGrid(:,3))];
                    
                    nPx = (diff(xLimsData) / cellSize.target) + 1;
                    nPz = (diff(zLimsData) / cellSize.target) + 1;
                    
                    sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
                    sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));
                    
                    cellSize.(planes{i}).x = sizeX(1); clear sizeX;
                    cellSize.(planes{i}).y = cellSize.target;
                    cellSize.(planes{i}).z = sizeZ(1); clear sizeZ;
                    cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).z;
        
                    xOrig = uData.(planes{i}).positionGrid(:,1);
                    zOrig = uData.(planes{i}).positionGrid(:,3);
                    
                    [x, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                                    linspace(zLimsData(1), zLimsData(2), nPz));
                    
                    uData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
                    uData.(planes{i}).positionGrid(:,2) = yLimsData;
                    uData.(planes{i}).positionGrid(:,[1,3]) = [x(:), z(:)];
                    
                    uInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).u.mean, 'linear', 'none');
                    vInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).v.mean, 'linear', 'none');
                    wInterp = scatteredInterpolant(xOrig, zOrig, uData.(planes{i}).w.mean, 'linear', 'none');
                    
                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,3));

                case 'XY'
                    xLimsData = [min(uData.(planes{i}).positionGrid(:,1)); ...
                                 max(uData.(planes{i}).positionGrid(:,1))];
                    yLimsData = [min(uData.(planes{i}).positionGrid(:,2)); ...
                                 max(uData.(planes{i}).positionGrid(:,2))];
                    zLimsData = uData.(planes{i}).positionGrid(1,3);
                    
                    nPx = (diff(xLimsData) / cellSize.target) + 1;
                    nPy = (diff(yLimsData) / cellSize.target) + 1;
                    
                    sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
                    sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
                    
                    cellSize.(planes{i}).x = sizeX(1); clear sizeX;
                    cellSize.(planes{i}).y = sizeY(1); clear sizeY;
                    cellSize.(planes{i}).z = cellSize.target;
                    cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).y;
        
                    xOrig = uData.(planes{i}).positionGrid(:,1);
                    yOrig = uData.(planes{i}).positionGrid(:,2);
                    
                    [x, y] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                                    linspace(yLimsData(1), yLimsData(2), nPy));
                    
                    uData.(planes{i}).positionGrid = zeros([height(x(:)),3]);
                    uData.(planes{i}).positionGrid(:,3) = zLimsData;
                    uData.(planes{i}).positionGrid(:,[1,2]) = [x(:), y(:)];
                    
                    uInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).u.mean, 'linear', 'none');
                    vInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).v.mean, 'linear', 'none');
                    wInterp = scatteredInterpolant(xOrig, yOrig, uData.(planes{i}).w.mean, 'linear', 'none');
                    
                    uData.(planes{i}).u.mean = uInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));
                    uData.(planes{i}).v.mean = vInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));
                    uData.(planes{i}).w.mean = wInterp(uData.(planes{i}).positionGrid(:,1), ...
                                                       uData.(planes{i}).positionGrid(:,2));

            end
            clear yOrig zOrig y z uInterp vInterp wInterp;

        end
        clear i;

end

disp(' ');

% Normalise Velocity
disp('    Normalising Velocity...');

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
    
    case 'C'
        
        if strcmp(campaignID, 'Varney')
            U = 40; % m/s
            
            for i = 1:height(planes)
                    uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
                    uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
                    uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
                    uData.(planes{i}).u.RMS = uData.(planes{i}).u.RMS / U;
                    uData.(planes{i}).v.RMS = uData.(planes{i}).v.RMS / U;
                    uData.(planes{i}).w.RMS = uData.(planes{i}).w.RMS / U;
            end
            
        end
        
end
        
%     case 'B'
%         
%         if contains(caseName, ["Run_Test", "Windsor"])
%             U = 40; % m/s
%             
%             for i = 1:height(planes)
%                 uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
%                 uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
%                 uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
%                 
%                 for j = 1:height(uData.(planes{i}).time)
%                     uData.(planes{i}).u{j} = uData.(planes{i}).u{j} / U;
%                     uData.(planes{i}).v{j} = uData.(planes{i}).v{j} / U;
%                     uData.(planes{i}).w{j} = uData.(planes{i}).w{j} / U;
%                     uData.(planes{i}).uPrime{j} = uData.(planes{i}).uPrime{j} / U;
%                     uData.(planes{i}).vPrime{j} = uData.(planes{i}).vPrime{j} / U;
%                     uData.(planes{i}).wPrime{j} = uData.(planes{i}).wPrime{j} / U;
%                 end
%                 
%             end
%             
%         end

disp(' ');

% Perform Blockage Correction
disp('    Performing Blockage Correction...');

switch format
    
    case 'A'
        
        if strcmp(campaignID, 'Windsor_Upstream_2023')
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
            
            for i = 1:height(planes)
                uData.(planes{i}).u.mean = uData.(planes{i}).u.mean * (At / (At - Am));
                uData.(planes{i}).v.mean = uData.(planes{i}).v.mean * (At / (At - Am));
                uData.(planes{i}).w.mean = uData.(planes{i}).w.mean * (At / (At - Am));
            end
            clear i;
            
        end

end

%     case 'B'
%         
%         if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
%             Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
%             At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
%             
%             for i = 1:height(planes)
%                 uData.(planes{i}).u.mean = uData.(planes{i}).u.mean * (At / (At - Am));
%                 uData.(planes{i}).v.mean = uData.(planes{i}).v.mean * (At / (At - Am));
%                 uData.(planes{i}).w.mean = uData.(planes{i}).w.mean * (At / (At - Am));
%                 
%                 for j = 1:height(uData.(planes{i}).time)
%                     uData.(planes{i}).u{j} = uData.(planes{i}).u{j} * (At / (At - Am));
%                     uData.(planes{i}).v{j} = uData.(planes{i}).v{j} * (At / (At - Am));
%                     uData.(planes{i}).w{j} = uData.(planes{i}).w{j} * (At / (At - Am));
%                     uData.(planes{i}).uPrime{j} = uData.(planes{i}).uPrime{j} * (At / (At - Am));
%                     uData.(planes{i}).vPrime{j} = uData.(planes{i}).vPrime{j} * (At / (At - Am));
%                     uData.(planes{i}).wPrime{j} = uData.(planes{i}).wPrime{j} * (At / (At - Am));
%                 end
%                 
%             end
%             
%         end

disp(' ');

% Calculate Vorticity
disp('    Calculating Vorticity...');

switch format

    case {'A', 'C'}

        for i = 1:height(planes)
        
            if contains(planes{i}, '_X_')
                orientation = 'YZ';
            elseif contains(planes{i}, '_Y_')
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
                    
                    [uData.(planes{i}).vorticity.mean, ~] = curl(y, z, v, w);
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean';
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean(:);
                    
                    % Omit Values Close to the Ground Plane
                    if strcmp(campaignID, 'Windsor_fullScale')
                        uData.(planes{i}).vorticity.mean(uData.(planes{i}).positionGrid(:,3) < 32e-3) = 0;
                    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                        uData.(planes{i}).vorticity.mean(uData.(planes{i}).positionGrid(:,3) < 8e-3) = 0;
                    end
                    
                case 'XZ'
                    gridShape = [height(unique(uData.(planes{i}).positionGrid(:,1))), ...
                                 height(unique(uData.(planes{i}).positionGrid(:,3)))];
                    
                    x = reshape(uData.(planes{i}).positionGrid(:,1), gridShape)';
                    z = reshape(uData.(planes{i}).positionGrid(:,3), gridShape)';
                    
                    u = reshape(uData.(planes{i}).u.mean, gridShape)';
                    w = reshape(uData.(planes{i}).w.mean, gridShape)';
                    
                    [uData.(planes{i}).vorticity.mean, ~] = curl(x, z, u, w);
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean';
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean(:);
                    
                    % Omit Values Close to the Ground Plane
                    if strcmp(campaignID, 'Windsor_fullScale')
                        uData.(planes{i}).vorticity.mean(uData.(planes{i}).positionGrid(:,3) < 32e-3) = 0;
                    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                        uData.(planes{i}).vorticity.mean(uData.(planes{i}).positionGrid(:,3) < 8e-3) = 0;
                    end
                    
                case 'XY'
                    gridShape = [height(unique(uData.(planes{i}).positionGrid(:,1))), ...
                                 height(unique(uData.(planes{i}).positionGrid(:,2)))];
                    
                    x = reshape(uData.(planes{i}).positionGrid(:,1), gridShape)';
                    y = reshape(uData.(planes{i}).positionGrid(:,2), gridShape)';
                    
                    u = reshape(uData.(planes{i}).u.mean, gridShape)';
                    v = reshape(uData.(planes{i}).v.mean, gridShape)';
                    
                    [uData.(planes{i}).vorticity.mean, ~] = curl(x, y, u, v);
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean';
                    uData.(planes{i}).vorticity.mean = uData.(planes{i}).vorticity.mean(:);
                    
            end
            
            % Omit Artificially High Values
            vortLims = prctile(uData.(planes{i}).vorticity.mean, [1, 99]);
            uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean < vortLims(1)) = vortLims(1);
            uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean > vortLims(2)) = vortLims(2);

            % Normalise
            uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean < 0) = ...
            rescale(uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean < 0), -1, 0);
            uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean > 0) = ...
            rescale(uData.(planes{i}).vorticity.mean(uData.(planes{i}).vorticity.mean > 0), 0, 1);
            uData.(planes{i}).vorticity.mean((uData.(planes{i}).vorticity.mean > -0.25) & ...
            (uData.(planes{i}).vorticity.mean < 0.25)) = 0;
            clear x y z u v w vortLims;
            
        end
        clear i;
        
end

disp(' ');

% Normalise Coordinate System
disp('    Normalising Spatial Dimensions...');

if normDims
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
    end
    clear i parts;
    
    xDims = xDims / normLength;
    yDims = yDims / normLength;
    zDims = zDims / normLength;
    
    cellSize.target = cellSize.target / normLength;
    
    for i = 1:height(planes)
        cellSize.(planes{i}).x = cellSize.(planes{i}).x / normLength;
        cellSize.(planes{i}).y = cellSize.(planes{i}).y / normLength;
        cellSize.(planes{i}).z = cellSize.(planes{i}).z / normLength;
        cellSize.(planes{i}).area = cellSize.(planes{i}).area / (normLength^2);

        uData.(planes{i}).positionGrid = uData.(planes{i}).positionGrid / normLength;
    end
    
end

% if normDims
% 
%     switch format
% 
%         case {'A', 'C'}
%             
%             for i = 1:height(planes)
%                 uData.(planes{i}).positionGrid = round((uData.(planes{i}).positionGrid / normLength), spacePrecision);
%                 
%                 [uData.(planes{i}).positionGrid, index] = unique(uData.(planes{i}).positionGrid, 'rows', 'stable');
% 
%                 uData.(planes{i}).u.mean = uData.(planes{i}).u.mean(index);
%                 uData.(planes{i}).v.mean = uData.(planes{i}).v.mean(index);
%                 uData.(planes{i}).w.mean = uData.(planes{i}).w.mean(index);
%             end
%             clear i;
%     
%     end
% 
% end

%%%%

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Map(s)? [y/n]: ', 's');

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

    case {'A', 'C'}
        % Select Plane(s) of Interest
        valid = false;
        while ~valid
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'multiple', ...
                                     'name', 'Select Variable(s) to Plot', ...
                                     'listString', planes);
        
            if ~valid
                disp(    'WARNING: No Planes Selected');
            end
        end
        clear valid;
        
        plotPlanes = planes(index);
        
        plotRMS = false;
        plotInst = false;
    
    case 'B'
        valid = false;
        while ~valid
            disp(' ');

            selection = input('Plot RMS Map? [y/n]: ', 's');

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

            selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

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

end

if plotMean || plotRMS || plotInst
    
    valid = false;
    while ~valid
        disp(' ');

        selection = input('Plot Vorticity Map(s)? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            plotVort = false;

            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            plotVort = true;

            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end

    end
    clear valid;

end

disp(' ');
disp(' ');


%% Present Velocity Data

disp('Map Presentation');
disp('-----------------');

disp(' ');

switch format

    case {'A', 'C'}
        
        if plotMean
            
            if normDims
                spatialRes = 0.5e-3 / normLength;
            else
                
                if strcmp(campaignID, 'Windsor_fullScale')
                    spatialRes = 2e-3;
                elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                    spatialRes = 0.5e-3;
                else
                    spatialRes = 0.5e-3;
                end

            end
            
            component = [];
            mapPerim = [];
            nPlanes = 1;
            planeNo = 1;
            figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
            
            for i = 1:height(plotPlanes)
                disp(['    Presenting ', plotPlanes{i}, '...']);
                
                if contains(plotPlanes{i}, '_X_')
                    orientation = 'YZ';
                elseif contains(plotPlanes{i}, '_Y_')
                    orientation = 'XZ';
                else
                    orientation = 'XY';
                end
        
                positionData = uData.(plotPlanes{i}).positionGrid;
                vectorData = [uData.(plotPlanes{i}).u.mean, uData.(plotPlanes{i}).v.mean, uData.(plotPlanes{i}).w.mean];
        
                switch orientation
        
                    case 'YZ'
                        xLimsPlot = [0.3; 4.6257662];
                        yLimsPlot = [-0.5; 0.5];
                        zLimsPlot = [0; 0.5];
        
                    case {'XZ', 'XY'}
                        xLimsPlot = [0.3; 1.2];
                        yLimsPlot = [-0.3; 0.3];
                        zLimsPlot = [0; 0.5];
                
                end
                
                if ~normDims

                    if strcmp(campaignID, 'Windsor_fullScale')
                        xLimsPlot = round((xLimsPlot * 4.176), spacePrecision);
                        yLimsPlot = round((yLimsPlot * 4.176), spacePrecision);
                        zLimsPlot = round((zLimsPlot * 4.176), spacePrecision);
                    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                        xLimsPlot = round((xLimsPlot * 1.044), spacePrecision);
                        yLimsPlot = round((yLimsPlot * 1.044), spacePrecision);
                        zLimsPlot = round((zLimsPlot * 1.044), spacePrecision);
                    end

                end
                
                switch orientation
        
                    case 'YZ'
                        xLimsData = uData.(plotPlanes{i}).positionGrid(1,1);
                        yLimsData = yLimsPlot;
                        zLimsData = zLimsPlot;
        
                    case 'XZ'
                        xLimsData = xLimsPlot;
                        yLimsData = uData.(plotPlanes{i}).positionGrid(1,2);
                        zLimsData = zLimsPlot;
        
                    case 'XY'
                        xLimsData = xLimsPlot;
                        yLimsData = yLimsPlot;
                        zLimsData = uData.(plotPlanes{i}).positionGrid(1,3);
        
                end

                switch orientation
                    
                    case 'YZ'
                        nComponents = 3;
                        
                    case {'XZ', 'XY'}
                        nComponents = 2;
                        
                end
                
                figName = [plotPlanes{i}, '_', caseID];
                cMap = viridis(32);
                streamlines = true;
                cLims = [0, 1];
                                    
                [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                       xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                                       mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                       streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                       zLimsPlot, normDims, figSave);
                
                if plotVort
                    scalarData = uData.(plotPlanes{i}).vorticity.mean;
                    figName = [plotPlanes{i}, '_Vorticity_', caseID];
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
            
            disp(' ');
        end

end

% switch format
% 
%     case 'B'
% 
%         if plotInst
%             for i = 1:height(plotPlanes)
%                 figHold = fig;
%                 
%                 for j = 1:nFrames
%                     
%                     if j ~= 1
%                         clf(fig)
%                         fig = figHold;
%                     end
% 
%                     vectorData = [uData.(plotPlanes{i}).u{j}, uData.(plotPlanes{i}).v{j}, uData.(plotPlanes{i}).w{j}];
%                     figTime = num2str(uData.(plotPlanes{i}).time(j), ['%.', num2str(timePrecision), 'f']);
%                     figName = [caseName, '_', plotPlanes{i}, '_T', erase(figTime, '.')];
%                     figSubtitle = [figTime, ' \it{s}'];
%                     
%                     [fig, planeNo] = plotPlanarVectorField_legacy(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
%                                                                   vectorData, nComponents, component, mapPerim, ...
%                                                                   nPlanes, planeNo, fig, figName, cMap, geometry, ...
%                                                                   streamlines, xDims, yDims, zDims, figTitle, figSubtitle, ...
%                                                                   cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise, figSave);
%                 end
%                 
%             end
%             
%         end
%         
% end


%% Local Functions

function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
