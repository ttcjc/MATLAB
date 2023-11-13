%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Velocity Processing v5.0

figSave = false; % Save .fig File(s)

normDims = true; % Normalise Spatial Dimensions

disp('===============================');
disp('Planar Velocity Processing v5.0');
disp('===============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Improve Versatility
% v3.0 - Moved Plotting to velocityPlots.m
% v4.0 - Rewrite to Support ParaView, Probe and Experimental Planar Data
% v5.0 - Update To Improve Consistency of Structures Across Repository


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
        [campaignID, caseID, uData] = initialisePVdata(saveLocation, 'U');
        
    case {'B', 'C'}
        error('NYI');
        
end

%     case 'B'
%         [caseName, dataID, velData, sampleInterval, timePrecision, geometry, ...
%          xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(saveLocation, 'planar', normalise, nProc);    

%     case 'C'
%         [caseName, velData, geometry, ...
%          xDims, yDims, zDims, spacePrecision] = initialiseExpData('U', normalise);

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(normDims);

disp(' ');
disp(' ');


%% Data Formatting

planes = fieldnames(uData);

% Adjust Data Origin
switch format

    case 'A'

        if strcmp(campaignID, 'Windsor_Upstream_2023')
            
            for i = 1:height(planes)
                uData.(planes{i}).positionGrid(:,1) = uData.(planes{i}).positionGrid(:,1) + 1.325;
            end

        end

end
        
%     case 'B'
%         
%         if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
%             
%             for i = 1:height(planes)
%                 
%                 if normalise
%                     velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);
%                     
%                     switch velData.(planes{i}).planeOrientation
%                         
%                         case 'YZ'
%                             velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + round((1.325 / 1.044), spacePrecision);
%                             velData.(planes{i}).xLims = velData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
%                             
%                         case {'XZ', 'XY'}
%                             velData.(planes{i}).xLims = velData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
%                             
%                     end
%                     
%                 else
%                     velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>
%                     
%                     switch velData.(planes{i}).planeOrientation
%                         
%                         case 'YZ'
%                             velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + 1.325;
%                             velData.(planes{i}).xLims = velData.(planes{i}).xLims + 1.325;
%                         
%                         case {'XZ', 'XY'}
%                             velData.(planes{i}).xLims = velData.(planes{i}).xLims + 1.325;
%                     
%                     end
%                     
%                 end
% 
%             end
%             
%         

% Normalise Coordinate System
if normDims

    switch format

        case 'A'
            
            for i = 1:height(planes)
                uData.(planes{i}).positionGrid = round((uData.(planes{i}).positionGrid / normLength), spacePrecision);
                
                [uData.(planes{i}).positionGrid, index] = unique(uData.(planes{i}).positionGrid, 'rows', 'stable');

                uData.(planes{i}).u.mean = uData.(planes{i}).u.mean(index);
                uData.(planes{i}).v.mean = uData.(planes{i}).v.mean(index);
                uData.(planes{i}).w.mean = uData.(planes{i}).w.mean(index);
            end
    
    end

end

% Normalise Velocity
switch format
    
    case 'A'
        
        if strcmp(campaignID, 'Windsor_fullScale')
            U = 22.22;
        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            U = 40;
        else
            U = 1;
        end
        
        for i = 1:height(planes)
            uData.(planes{i}).u.mean = uData.(planes{i}).u.mean / U;
            uData.(planes{i}).v.mean = uData.(planes{i}).v.mean / U;
            uData.(planes{i}).w.mean = uData.(planes{i}).w.mean / U;
        end

end
        
%     case 'B'
%         
%         if contains(caseName, ["Run_Test", "Windsor"])
%             U = 40; % m/s
%             
%             for i = 1:height(planes)
%                 velData.(planes{i}).u.mean = velData.(planes{i}).u.mean / U;
%                 velData.(planes{i}).v.mean = velData.(planes{i}).v.mean / U;
%                 velData.(planes{i}).w.mean = velData.(planes{i}).w.mean / U;
%                 
%                 for j = 1:height(velData.(planes{i}).time)
%                     velData.(planes{i}).u{j} = velData.(planes{i}).u{j} / U;
%                     velData.(planes{i}).v{j} = velData.(planes{i}).v{j} / U;
%                     velData.(planes{i}).w{j} = velData.(planes{i}).w{j} / U;
%                     velData.(planes{i}).uPrime{j} = velData.(planes{i}).uPrime{j} / U;
%                     velData.(planes{i}).vPrime{j} = velData.(planes{i}).vPrime{j} / U;
%                     velData.(planes{i}).wPrime{j} = velData.(planes{i}).wPrime{j} / U;
%                 end
%                 
%             end
%             
%         end
        
%     case 'C'
%         
%         if strcmp(caseName, 'Varney')
%             U = 40; % m/s
%             
%             for i = 1:height(planes)
%                     velData.(planes{i}).u.mean = velData.(planes{i}).u.mean / U;
%                     velData.(planes{i}).v.mean = velData.(planes{i}).v.mean / U;
%                     velData.(planes{i}).w.mean = velData.(planes{i}).w.mean / U;
%                     velData.(planes{i}).uRMS = velData.(planes{i}).uRMS / U;
%                     velData.(planes{i}).vRMS = velData.(planes{i}).vRMS / U;
%                     velData.(planes{i}).wRMS = velData.(planes{i}).wRMS / U;
%             end
%             
%         end

% Perform Blockage Correction
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

% Map Raw Data Onto Adjusted Grid
if normDims
    cellSize.target = 1e-3;
else

    if strcmp(campaignID, 'Windsor_fullScale')
        cellSize.target = 4e-3;
    else
        cellSize.target = 1e-3;
    end

end

switch format

    case 'A'

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
                    yLimsData = [min(uData.(planes{i}).positionGrid(:,2)); max(uData.(planes{i}).positionGrid(:,2))];
                    zLimsData = [min(uData.(planes{i}).positionGrid(:,3)); max(uData.(planes{i}).positionGrid(:,3))];
        
                    cellSize.(planes{i}).x = cellSize.target;
                    cellSize.(planes{i}).y = (yLimsData(2) - yLimsData(1)) / ...
                             round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
                    cellSize.(planes{i}).z = (zLimsData(2) - zLimsData(1)) / ...
                             round(((zLimsData(2) - zLimsData(1)) / cellSize.target));
        
                    cellSize.(planes{i}).area = cellSize.(planes{i}).y * cellSize.(planes{i}).z;
        
                    yOrig = uData.(planes{i}).positionGrid(:,2);
                    zOrig = uData.(planes{i}).positionGrid(:,3);
                    
                    [y, z] = ndgrid(yLimsData(1):cellSize.(planes{i}).y:yLimsData(2), ...
                                    zLimsData(1):cellSize.(planes{i}).z:zLimsData(2));
                    
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
                    xLimsData = [min(uData.(planes{i}).positionGrid(:,1)); max(uData.(planes{i}).positionGrid(:,1))];
                    yLimsData = uData.(planes{i}).positionGrid(1,2);
                    zLimsData = [min(uData.(planes{i}).positionGrid(:,3)); max(uData.(planes{i}).positionGrid(:,3))];
        
                    cellSize.(planes{i}).x = (xLimsData(2) - xLimsData(1)) / ...
                             round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
                    cellSize.(planes{i}).y = cellSize.target;
                    cellSize.(planes{i}).z = (zLimsData(2) - zLimsData(1)) / ...
                             round(((zLimsData(2) - zLimsData(1)) / cellSize.target));
        
                    cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).z;
        
                    xOrig = uData.(planes{i}).positionGrid(:,1);
                    zOrig = uData.(planes{i}).positionGrid(:,3);
                    
                    [x, z] = ndgrid(xLimsData(1):cellSize.(planes{i}).x:xLimsData(2), ...
                                    zLimsData(1):cellSize.(planes{i}).z:zLimsData(2));
                    
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
                    xLimsData = [min(uData.(planes{i}).positionGrid(:,1)); max(uData.(planes{i}).positionGrid(:,1))];
                    yLimsData = [min(uData.(planes{i}).positionGrid(:,2)); max(uData.(planes{i}).positionGrid(:,2))];
                    zLimsData = uData.(planes{i}).positionGrid(1,3);

                    cellSize.(planes{i}).x = (xLimsData(2) - xLimsData(1)) / ...
                             round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
                    cellSize.(planes{i}).y = (yLimsData(2) - yLimsData(1)) / ...
                             round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
                    cellSize.(planes{i}).z = cellSize.target;
        
                    cellSize.(planes{i}).area = cellSize.(planes{i}).x * cellSize.(planes{i}).y;
        
                    xOrig = uData.(planes{i}).positionGrid(:,1);
                    yOrig = uData.(planes{i}).positionGrid(:,2);
                    
                    [x, y] = ndgrid(xLimsData(1):cellSize.(planes{i}).x:xLimsData(2), ...
                                    yLimsData(1):cellSize.(planes{i}).y:yLimsData(2));
                    
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

end


%% Vorticity Calculations



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

    case 'A'
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
        
    otherwise
        plotRMS = false;
        plotInst = false;
        
end

disp(' ');
disp(' ');


%% Present Velocity Data

disp('Map Presentation');
disp('-----------------');

disp(' ');

switch format

    case 'A'
        
        if plotMean

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
                
                if normDims
                    spatialRes = 0.5e-3;
                else
                    
                    if strcmp(campaignID, 'Windsor_fullScale')
                        spatialRes = 2e-3;
                    else
                        spatialRes = 0.5e-3;
                    end
                    
                end
        
                switch orientation
        
                    case 'YZ'
                        
                        if normDims
                            xLimsPlot = [0.3; 5.5];
                            yLimsPlot = [-0.4; 0.4];
                            zLimsPlot = [0; 0.6];
                        else
                            
                            if strcmp(campaignID, 'Windsor_fullScale')
                                xLimsPlot = [1.2528; 22.968];
                                yLimsPlot = [-2.088 2.088];
                                zLimsPlot = [0; 2.5056];
                            else
                                xLimsPlot = [0.3132; 5.742];
                                yLimsPlot = [-0.522; 0.522];
                                zLimsPlot = [0; 0.6264];
                            end
                        
                        end
        
                    case {'XZ', 'XY'}
                
                        if normDims
                            xLimsPlot = [0.3; 1.3];
                            yLimsPlot = [-0.4; 0.4];
                            zLimsPlot = [0; 0.6];
                        else
                            
                            if strcmp(campaignID, 'Windsor_fullScale')
                                xLimsPlot = [1.2528; 5.4288];
                                yLimsPlot = [-2.088 2.088];
                                zLimsPlot = [0; 2.5056];
                            else
                                xLimsPlot = [0.3132; 1.3572];
                                yLimsPlot = [-0.522; 0.522];
                                zLimsPlot = [0; 0.6264];
                            end
                        
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

                nComponents = 3;
                component = [];
                mapPerim = [];
                nPlanes = 1;
                planeNo = 1;
                figName = [plotPlanes{i}, '_', caseID];
                cMap = viridis(32);
                streamlines = true;
                figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
                cLims = [0, 1];
                                    
                [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                       xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                                       mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                       streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                       zLimsPlot, normDims, figSave);
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