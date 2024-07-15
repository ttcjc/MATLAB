run preamble;

%#ok<*UNRCH>

normDensity = false;

normDims = true;

figSave = false;

plotMean = true;

plotRMS = true;

plotInst = false;

format = 'B';

startFrame = 1;

endFrame = 5;


%%

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData(geoLoc);

if normDims

    parts = fieldnames(geometry);
    for i = 1:height(parts)
        geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
    end
    clear i parts;
    
end


%%

load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P6_109/T1002_T3200_F50_D20_D400_cumulative.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P10_285/T1002_T3200_F50_D20_D400_cumulative.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P14_461/T1002_T3200_F50_D20_D400_cumulative.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P18_637/T1002_T3200_F50_D20_D400_cumulative.mat');

nTimes = height(mapData.time);


%%

mapDataVars = fieldnames(mapData);
nonFieldVars = {'positionGrid'; 'time'; 'CoM'};
fieldVars = setdiff(mapDataVars, nonFieldVars);
clear mapDataVars nonFieldVars;

valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'multiple', ...
                             'name', 'Select Variable(s) to Plot', ...
                             'listString', fieldVars);

    if ~valid
        disp('    WARNING: No Mapping Variable Selected');
    end

end
clear valid;

plotVars = fieldVars(index); clear fieldVars;


%%

if normDims
    cellSize.target = cellSize.target / normLength;
    cellSize.y = cellSize.y / normLength;
    cellSize.z = cellSize.z / normLength;
    cellSize.area = cellSize.area / (normLength^2);

    mapData.positionGrid = mapData.positionGrid / normLength;
end

if strcmp(campaignID, 'Windsor_Upstream_2023')
    refValue = 8.996259860381801e-05;
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
    refValue = 0.137170217583008;
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
    refValue = 0.710023688287996;
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
    refValue = 0.226673182043117;
elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
    refValue = 0.441874462356553;
else
    refValue = prctile(full(mapData.areaDensity.mean(mapData.areaDensity.mean > 0)), 99);
end

if normDensity
    
    mapData.areaDensity.mean = mapData.areaDensity.mean / refValue;
    mapData.areaDensity.RMS = mapData.areaDensity.RMS / refValue;
    
    for i = 1:nTimes
        mapData.areaDensity.inst{i} = mapData.areaDensity.inst{i} / refValue;
        mapData.areaDensity.prime{i} = mapData.areaDensity.prime{i} / refValue;
    end
    clear i;
    
end


%%

clc;
close all;

mapPerim = [];

xLimsData = mapData.positionGrid(1,1) * normLength;
yLimsData = [-0.5; 0.5] * normLength;
zLimsData = [0; 0.5] * normLength;

if normDims
    xLimsData = xLimsData / normLength;
    yLimsData = yLimsData / normLength;
    zLimsData = zLimsData / normLength;
end

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    positionData = mapData.positionGrid;
    
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
    
    % Offset Particles From Surface to Improve Visibility
    switch format
        
        case 'A'
            xLimsData = xLimsData + spatialRes;
            positionData(:,1) = xLimsData;
            
    end
    
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    
    switch format

        case 'A'
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.25; 0.25];
            zLimsPlot = [0; 0.4];
            
        case 'B'
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];
            
    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end

if plotMean
    
    for i = 1:height(plotVars)        
        scalarData = full(mapData.(plotVars{i}).mean);
        
        switch format
            
            case 'A'
                figName = ['Average_Base_', plotVars{i}, '_', caseID];
                
            case 'B'
                figName = ['Average_', planeID, '_', plotVars{i}, '_', caseID];
                
        end
        
        switch format
            
            case 'A'
                contourlines = [];
                
            case 'B'
                
                if strcmp(plotVars{i}, 'areaDensity')
                    
                    if normDensity
                        contourlines = [0.02; 0.02];
                    else
                        contourlines = [0.02; 0.02] * refValue;
                    end
                    
                end
                
                if ~exist('contourlines', 'var')
                    contourlines = [];
                end
                
        end
        
        figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
            
            if strcmp(campaignID, 'Windsor_fullScale')
                cLims = [20; 210];
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                cLims = [0; 40];
            end
            
        elseif strcmp(plotVars{i}, 'areaDensity')
                
            switch format

                case 'A'

                    if normDensity
                        cLims = [0; 1.05];
                    else
                        %
                    end

                case 'B'

                    if normDensity
                        cLims = [0; 1.05];
                    else
                        
                        if strcmp(campaignID, 'Windsor_fullScale')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 885e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                cLims = [0; 215e-3];
                            elseif strcmp(planeID, 'X_P14_461')
                                cLims = [0; 115e-3];
                            elseif strcmp(planeID, 'X_P18_637')
%                                 cLims = [0; 55e-3];
                                cLims = [0; 90e-3];
                            end
                            
                        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                            %
                        end     
                        
                    end

            end
            
        end
            
        if ~exist('cLims', 'var')
            cLims = [0; max(scalarData)];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
        
        clear cLims;
    end
    clear i;
    
end


if plotRMS
    
    for i = 1:height(plotVars)        
        scalarData = full(mapData.(plotVars{i}).RMS);
        
        switch format
            
            case 'A'
                figName = ['RMS_Base_', plotVars{i}, '_', caseID];
                
            case 'B'
                figName = ['RMS_', planeID, '_', plotVars{i}, '_', caseID];
                
        end
        
        contourlines = [];
        figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
        
        if strcmp(plotVars{i}, 'areaDensity')
                
            switch format

                case 'A'

                    if normDensity
                        %
                    else
                        %
                    end

                case 'B'

                    if normDensity
                        %
                    else
                        
                        if strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 110e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 25e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 640e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 50e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 210e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 25e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 380e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 75e-3];
                            end
                            
                        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                            %
                        end     
                        
                    end

            end
            
        end
            
        if ~exist('cLims', 'var')
            cLims = [0; max(scalarData)];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
        
        clear cLims;
    end
    clear i;
    
end

if plotInst
    
    for i = 1:height(plotVars)        
        contourlines = [];
        
        if any(strcmp(plotVars{i}, {'d10', 'd20', 'd30', 'd32'}))
            
            if strcmp(campaignID, 'Windsor_fullScale')
                cLims = [20; 400];
            else
                cLims = [0; 150];
            end
            
        elseif strcmp(plotVars{i}, 'areaDensity')
            
            switch format

                case 'A'

                    if normDensity
                        %
                    else
                        %
                    end

                case 'B'

                    if normDensity
                        %
                    else
                        
                        if strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 780e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 210e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 3590e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 385e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 1375e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 210e-3];
                            end
                            
                        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase20deg')
                            
                            if strcmp(planeID, 'X_P6_109')
                                cLims = [0; 2380e-3];
                            elseif strcmp(planeID, 'X_P10_285')
                                %
                            elseif strcmp(planeID, 'X_P14_461')
                                %
                            elseif strcmp(planeID, 'X_P18_637')
                                cLims = [0; 595e-3];
                            end
                            
                        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                            %
                        end     
                        
                    end

            end
            
        end
        
        if ~exist('cLims', 'var') || isempty(cLims)
            instMax = zeros([nTimes,1]);

            for j = 1:nTimes
%                 instMax(j) = prctile(mapData.areaDensity.inst{j}, 99);
                instMax(j) = max(mapData.areaDensity.inst{j});
            end
            clear j;

%             cLims = [0; max(instMax)];
            cLims = [0; prctile(instMax, 99)];
        end
        
        figHold = fig;
        
        for j = startFrame:endFrame
            
            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            scalarData = full(mapData.(plotVars{i}).inst{j});
            figTime = num2str(mapData.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Inst_Base_', plotVars{i}, '_T', erase(figTime, '.'), '_', caseID];
                
                case 'B'
                    figName = ['Inst_', planeID, '_', plotVars{i}, '_T', erase(figTime, '.'), '_', caseID];
            end
            
            figTitle = ['{', figTime, ' \it{s}}'];
        
            [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                   planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                   refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                   zLimsPlot, normDims, figSave);
        end
        clear j;
        
        clear cLims;        
    end
    clear i;
    
end