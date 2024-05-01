run preamble;

%#ok<*UNRCH>

normDensity = true;

normDims = true;

format = 'C';

figSave = false;

plotMean = false;

plotInst = true;

startFrame = 100;

endFrame = 100;


%%

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData(geoLoc);


%%

% load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/midWake/T12525_T40000_F400_D1_D147.mat');
load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/farWake/T1002_T3200_F50_D20_D400.mat');

nTimes = height(volumeData.time);

%%

switch format
    
    case 'A' % 1 L
        xLimsData = [-0.537116858237548; 1.462883141762452] * normLength;
        yLimsData = [-0.4; 0.4] * normLength;
        zLimsData = [0; 0.4] * normLength;
        
    case 'B' % 2 L
        xLimsData = [-0.537116858237548; 2.462883141762452] * normLength;
        yLimsData = [-0.5; 0.5] * normLength;
        zLimsData = [0; 0.5] * normLength;
        
    case 'C' % 4 L
        xLimsData = [-0.537116858237548; 4.462883141762452] * normLength;
        yLimsData = [-0.6; 0.6] * normLength;
        zLimsData = [0; 0.6] * normLength;

end

plane1L = [
           (xDims(2) + (1 * normLength)), yLimsData(1), zLimsData(1);
           (xDims(2) + (1 * normLength)), yLimsData(1), zLimsData(2);
           (xDims(2) + (1 * normLength)), yLimsData(2), zLimsData(2);
           (xDims(2) + (1 * normLength)), yLimsData(2), zLimsData(1);
           (xDims(2) + (1 * normLength)), yLimsData(1), zLimsData(1);
          ];

plane2L = [
           (xDims(2) + (2 * normLength)), yLimsData(1), zLimsData(1);
           (xDims(2) + (2 * normLength)), yLimsData(1), zLimsData(2);
           (xDims(2) + (2 * normLength)), yLimsData(2), zLimsData(2);
           (xDims(2) + (2 * normLength)), yLimsData(2), zLimsData(1);
           (xDims(2) + (2 * normLength)), yLimsData(1), zLimsData(1);
          ];

plane3L = [
           (xDims(2) + (3 * normLength)), yLimsData(1), zLimsData(1);
           (xDims(2) + (3 * normLength)), yLimsData(1), zLimsData(2);
           (xDims(2) + (3 * normLength)), yLimsData(2), zLimsData(2);
           (xDims(2) + (3 * normLength)), yLimsData(2), zLimsData(1);
           (xDims(2) + (3 * normLength)), yLimsData(1), zLimsData(1);
          ];

plane4L = [
           (xDims(2) + (4 * normLength)), yLimsData(1), zLimsData(1);
           (xDims(2) + (4 * normLength)), yLimsData(1), zLimsData(2);
           (xDims(2) + (4 * normLength)), yLimsData(2), zLimsData(2);
           (xDims(2) + (4 * normLength)), yLimsData(2), zLimsData(1);
           (xDims(2) + (4 * normLength)), yLimsData(1), zLimsData(1);
          ];


%%

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
    cellSize.x = cellSize.x / normLength;
    cellSize.y = cellSize.y / normLength;
    cellSize.z = cellSize.z / normLength;
    cellSize.volume = cellSize.volume / (normLength^3);

    volumeData.positionGrid = volumeData.positionGrid / normLength;
    
    xLimsData = xLimsData / normLength;
    yLimsData = yLimsData / normLength;
    zLimsData = zLimsData / normLength;
    
    plane1L = plane1L / normLength;
    plane2L = plane2L / normLength;
    plane3L = plane3L / normLength;
    plane4L = plane4L / normLength;

end

if normDensity

    if strcmp(campaignID, 'Windsor_Upstream_2023')
        refValue = 0.002551862743441; % Windsor_SB_wW_Upstream_SC
    else
        refValue = full(volumeData.density.mean); refValue = prctile(refValue(refValue > 0), 99);
    end

    volumeData.density.mean = volumeData.density.mean / refValue;

    for i = 1:nTimes
        volumeData.density.inst{i} = volumeData.density.inst{i} / refValue;
    end
    clear i;
    
end


%%

clc;
close all;

if plotMean || plotInst
    gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
                 height(unique(volumeData.positionGrid(:,2))), ...
                 height(unique(volumeData.positionGrid(:,3)))];
             
    spatialRes = cellSize.target / 2;
    xOrig = reshape(volumeData.positionGrid(:,1), gridShape);
    yOrig = reshape(volumeData.positionGrid(:,2), gridShape);
    zOrig = reshape(volumeData.positionGrid(:,3), gridShape);
    POD = false;
    nSurfaces = 1;
    surfaceNo = 1;
    
    if strcmp(campaignID, 'Windsor_fullScale')
        
        if strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
            cMap = graphColours(4);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
            cMap = graphColours(5);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
            cMap = graphColours(6);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
            cMap = graphColours(7);
        end
        
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        
        if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
            cMap = graphColours(1);
        elseif strcmp(caseID, 'Windsor_ST_wW_Upstream_SC')
            cMap = graphColours(2);
        elseif strcmp(caseID, 'Windsor_RSST_wW_Upstream_SC')
            cMap = graphColours(3);
        end
        
    end
    
    if ~exist('cMap', 'var')
        cMap = graphColours(1);
    end
    
    viewAngle = [30, 30];
    
    switch format

        case 'A' % 1 L
            xLimsPlot = [-0.637116858237548; 1.562883141762452];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];

        case 'B' % 2 L
            xLimsPlot = [-0.637116858237548; 2.562883141762452];
            yLimsPlot = [-0.6; 0.6];
            zLimsPlot = [0; 0.6];

        case 'C' % 4 L
            xLimsPlot = [-0.637116858237548; 4.562883141762452];
            yLimsPlot = [-0.7; 0.7];
            zLimsPlot = [0; 0.7];

    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end
    
if plotMean
    
    if normDensity
        fieldData = reshape((full(volumeData.density.mean) * 100), gridShape);
        isoValue = [25; 2];
    else
        fieldData = reshape(full(volumeData.density.mean), gridShape);
        isoValue = [0.25; 0.02] * max(full(volumeData.density.mean));
    end
    
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    multiView = false;
    
    for i = 1:height(isoValue)
        
        switch format

            case 'A'
                figName = ['NW_Average_Density_', num2str(isoValue(i)), '%_', caseID];

            case 'B'
                figName = ['MW_Average_Density_', num2str(isoValue(i)), '%_', caseID];

            case 'C'
                figName = ['FW_Average_Density_', num2str(isoValue(i)), '%_', caseID];

        end
        
        [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                           xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                           fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                           multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
        
    end
    clear i;
    
end

if plotInst
    
    if normDensity
        isoValue = [50; 10];
    else
        isoValue = 0.5 * max(full(volumeData.densivty.mean));
    end
    
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            if normDensity
                fieldData = reshape((full(volumeData.density.inst{j}) * 100), gridShape);
            else
                fieldData = reshape(full(volumeData.density.inst{j}), gridShape);
            end
            
            figTime = num2str(volumeData.time(j), ['%.', num2str(timePrecision), 'f']);

            switch format

                case 'A'
                    figName = ['NW_Inst_Density_', num2str(isoValue(i)), '%_T'...
                               erase(figTime, '.'), '_', caseID];
                    
                case 'B'
                    figName = ['MW_Inst_Density_', num2str(isoValue(i)), '%_T'...
                               erase(figTime, '.'), '_', caseID];
                    
                case 'C'
                    figName = ['FW_Inst_Density_', num2str(isoValue(i)), '%_T'...
                               erase(figTime, '.'), '_', caseID];

            end

            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                               xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                               fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                               multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
            
            hold on;
            
            plot3(plane1L(:,1), plane1L(:,2), plane1L(:,3), 'lineStyle', '-', ...
                                                            'lineWidth', 2, ...
                                                            'color', graphColours(2));
            plot3(plane2L(:,1), plane2L(:,2), plane2L(:,3), 'lineStyle', '-', ...
                                                            'lineWidth', 2, ...
                                                            'color', graphColours(2));
            plot3(plane3L(:,1), plane3L(:,2), plane3L(:,3), 'lineStyle', '-', ...
                                                            'lineWidth', 2, ...
                                                            'color', graphColours(2));
            plot3(plane4L(:,1), plane4L(:,2), plane4L(:,3), 'lineStyle', '-', ...
                                                            'lineWidth', 2, ...
                                                            'color', graphColours(2));
                                                        
            print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
        end
        clear j;
        
    end
    clear i;
    
end