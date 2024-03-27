clc;
close all;

planePosition = 1.949 / normLength;

% yRay = yLimsData(1):(diff(yLimsData) / 5):yLimsData(2);
% zRay = zLimsData(1):(diff(zLimsData) / 5):zLimsData(2);
% 
% [yRay, zRay] = ndgrid(yRay, zRay);
% 
% positionsRay = [(ones(height(yRay(:)),1) * planePosition) yRay(:), zRay(:)]; clear yRay zRay;

% positionsRay = dsearchn(volumeData.positionGrid(:,[2,3]), ([0, 0.76] / normLength));
% positionsRay = [xDims(2), volumeData.positionGrid(positionsRay, [2,3])];

positionsPlane = [
                  planePosition, yLimsData(1), zLimsData(1);
                  planePosition, yLimsData(1), zLimsData(2);
                  planePosition, yLimsData(2), zLimsData(2);
                  planePosition, yLimsData(2), zLimsData(1);
                  planePosition, yLimsData(1), zLimsData(1);
                 ];

positionsSensor = [18.637, -1.945, 0.76] / normLength;

if plotInst
    disp('Presenting Instantaneous Volume Field...');
    
    isoValue = 0.5;
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            if normDensity
                fieldData = reshape(full(volumeData.density.inst{j}), gridShape);
            else
                fieldData = reshape((full(volumeData.density.inst{j}) / ...
                                    prctile(full(volumeData.density.mean), 99)), gridShape);
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
            
            % Plot Ray Projection
            plot3(positionsSensor(:,1), positionsSensor(:,2), positionsSensor(:,3), 'lineStyle', 'none', ...
                                                                                    'marker', 'o', ...
                                                                                    'markerSize', 10, ...
                                                                                    'markerFaceColor', graphColours(2), ...
                                                                                    'markerEdgeColor', graphColours(2));
            
            plot3(positionsPlane(:,1), positionsPlane(:,2), positionsPlane(:,3), 'lineStyle', '-', ...
                                                                                 'lineWidth', 2, ...
                                                                                 'color', graphColours(2));
            
%             for k = 1:height(positionsRay)
% %                 plot3([positionsSensor(:,1); positionsRay(k,1)], ...
% %                       [positionsSensor(:,2); positionsRay(k,2)], ... 
% %                       [positionsSensor(:,3); positionsRay(k,3)], 'lineStyle', '-', ...
% %                                                                  'lineWidth', 2, ...
% %                                                                  'color', ([74, 24, 99, 127.5] / 255));
%                 
%                 plot3([positionsSensor(:,1); positionsRay(k,1)], ...
%                       [positionsSensor(:,2); positionsRay(k,2)], ... 
%                       [positionsSensor(:,3); positionsRay(k,3)], 'lineStyle', '-', ...
%                                                                  'lineWidth', 4, ...
%                                                                  'color', ([74, 24, 99] / 255));
%             end
%             clear k;
            
            for k = 1:4
                plot3([positionsSensor(:,1); positionsPlane(k,1)], ...
                      [positionsSensor(:,2); positionsPlane(k,2)], ... 
                      [positionsSensor(:,3); positionsPlane(k,3)], 'lineStyle', '-', ...
                                                                   'lineWidth', 2, ...
                                                                   'color', ([74, 24, 99, 127.5] / 255));
            end
            clear k;
            
%             print(gcf, [userpath, '/Output/Figures/Ray_Projection_Principle.png'], '-dpng', '-r300');
            print(gcf, [userpath, '/Output/Figures/LoS_Context.png'], '-dpng', '-r300');
%             print(gcf, [userpath, '/Output/Figures/Ray_Context.png'], '-dpng', '-r300'); 
        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end