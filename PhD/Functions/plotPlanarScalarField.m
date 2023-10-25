%% Planar Scalar Field Plotter v2.2
% ----
% Plots Previously Processed Planar Scalar Fields
% ----
% Usage: [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
%                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
%                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
%                                               refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
%                                               yLimsPlot, zLimsPlot, normalise, figSave);
% 
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        'positionData' -> Cartesian Positions of Data Points
%        'scalarData'   -> Field Value @ 'positionData' Points
%        'spatialRes'   -> Target Grid Spacing [Dimensions of 'positionData']
%        '*LimsData'    -> Contour Plot Limits [Dimensions of 'positionData']
%        'mapPerim'     -> Map Perimeter Used When Plotting a Plane of Arbitrary Shape [Dimensions of 'positionData']
%        'nPlanes'      -> Number of Planes in a Multi-Plane Figure
%        'planeNo'      -> Current Plane Number
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'cMap'         -> Colour Map
%        'geometry'     -> STL(s) to Include in Plot
%        'contourlines' -> Contour Lines at the Corresponding Values
%        'refPoint'     -> Reference Point (CoM, CoP etc.) [Dimensions of 'positionData']
%        'figTitle'     -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle'  -> Figure Title
%        'cLims'        -> Colour Map Limits
%        '*LimsPlot'    -> 3D Axes Limits [Dimensions of 'positionData']
%        'normalise'    -> Normalise Dimensions [True/False]
%        'figSave'     -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit (Expanded Functionality of 'contaminantPlots')
% v2.0 - Shifted to Using 'griddedInterpolant' and Added Support for Multi-Plane Figures)
% v2.1 - Rename and Minor Formatting Updates
% v2.1 - Added Spatial Resolution as an Input Variable
% v2.2 - Update To Ensure Consistent Figure Sizes


%% Main Function

function [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                                yLimsPlot, zLimsPlot, normalise, figSave)
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,2))), ...
                         height(unique(positionData(:,3)))];

            y = reshape(positionData(:,2), gridShape);
            z = reshape(positionData(:,3), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(y, z, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = spatialRes;
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / spatialRes));
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / spatialRes));
            
            [x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
            
            % Map Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(2,:,:) = interp(y(2,:,:), z(2,:,:));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(y, z, mapPerim(:,2), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,3)))];

            x = reshape(positionData(:,1), gridShape);
            z = reshape(positionData(:,3), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(x, z, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / spatialRes));
            cellSizeY = spatialRes;
            cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / spatialRes));

            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               (yLimsData - cellSizeY):cellSizeY:(yLimsData + cellSizeY), ...
                               zLimsData(1):cellSizeZ:zLimsData(2));
    
            % Map Scalar Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(:,2,:) = interp(x(:,2,:), z(:,2,:));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, z, mapPerim(:,1), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
        case 'XY'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,2)))];

            x = reshape(positionData(:,1), gridShape);
            y = reshape(positionData(:,2), gridShape);
            
            scalar = reshape(scalarData, gridShape);
            
            % Perform Interpolation
            interp = griddedInterpolant(x, y, scalar, 'linear', 'none');
            
            % Generate 3D Surface
            cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / spatialRes));
            cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / spatialRes));
            cellSizeZ = spatialRes;
            
            [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                               yLimsData(1):cellSizeY:yLimsData(2), ...
                               (zLimsData - cellSizeZ):cellSizeZ:(zLimsData + cellSizeZ));

            % Map Scalar Data on to 3D Surface
            scalar = zeros(size(x));
            scalar(:,:,2) = interp(x(:,:,2), y(:,:,2));
            
            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, y, mapPerim(:,1), mapPerim(:,2));
                indexMap = double(or(indexIn, indexOn));
                indexMap(indexMap == 0) = nan;
                
                scalar = scalar .* indexMap;
            end
            
    end
    
    % Initialise Figure
    if planeNo == 1                
        fig = fig + 1;
        set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                         'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
        pause(0.5);
        hold on;
        set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
                 'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
        lighting gouraud;
        colormap(cMap);

        % Plot Geometry
        if ~isempty(geometry)
            parts = fieldnames(geometry);

            for i = 1:height(parts)
                patch('faces', geometry.(parts{i}).faces, ...
                      'vertices', geometry.(parts{i}).vertices, ...
                      'faceColor', [0.5, 0.5, 0.5], ...
                      'edgeColor', [0.5, 0.5, 0.5], ...
                      'lineStyle', 'none');
            end
            clear i;

        end

    end
    
    % Present Data
    switch orientation
        
        case 'YZ'
            % Plot Scalar Field
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(scalar(2,:,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, xLimsData, [], [], contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(90, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([90, 0]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = (yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2));
                    yticks(tickData(2:(end-1)));
                    tickData = (zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2));
                    zticks(tickData(2:(end-1)));
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normalise
                        ylabel('{$y_{\ell}$}', 'interpreter', 'latex')
                        zlabel('{$z_{\ell}$}', 'interpreter', 'latex');
                    else
                        ylabel('{$y$ ($m$)}', 'interpreter', 'latex');
                        zlabel('{$z$ ($m$)}', 'interpreter', 'latex');
                    end
                
                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end
                
                tightInset = get(gca, 'TightInset');
                set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                           (tightInset(2) + 0.00625), ...
                                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
                pause(0.5);
                hold off;
                
                % Save Figure
                print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end
            
        case 'XZ'            
            % Plot Scalar Field
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(scalar(:,2,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, [], yLimsData, [], contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off
                    caxis(cLims);
                    view([0, 0]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:5));
                    tickData = [];
                    yticks(tickData);
                    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:5));
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normalise
                        xlabel('{$x_{\ell}$}', 'interpreter', 'latex')
                        zlabel('{$z_{\ell}$}', 'interpreter', 'latex');
                    else
                        xlabel('{$x$ ($m$)}', 'interpreter', 'latex');
                        zlabel('{$z$ ($m$)}', 'interpreter', 'latex');
                    end

                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end
                
                tightInset = get(gca, 'TightInset');
                set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                           (tightInset(2) + 0.00625), ...
                                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
                pause(0.5);
                hold off;
                
                % Save Figure
                print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end
            
        case 'XY'
            % Plot Scalar Field
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(scalar(:,:,2)), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);
            
            % Plot Contour Lines
            if ~isempty(contourlines)
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                scalar = permute(scalar, [2,1,3]);
                
                contours = contourslice(x, y, z, scalar, [], [], zLimsData, contourlines);
                set(contours, 'edgeColor', 'w', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Plot Reference Point
            if ~isempty(refPoint)
                scatter3(refPoint(1), refPoint(2), refPoint(3), 200, 'w', 'x', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([0, 90]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:5));
                    tickData = round((yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2)), 2);
                    yticks(tickData(2:5));
                    tickData = [];
                    zticks(tickData);
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normalise
                        xlabel('{$x_{\ell}$}', 'interpreter', 'latex')
                        ylabel('{$y_{\ell}$}', 'interpreter', 'latex');
                    else
                        xlabel('{$x$ ($m$)}', 'interpreter', 'latex');
                        ylabel('{$y$ ($m$)}', 'interpreter', 'latex');
                    end

                else
                    title(figTitle, 'color', ([254, 254, 254] / 255));
                    subtitle(figSubtitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    view([30, 30]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = [];
                    xticks(tickData);
                    tickData = [];
                    yticks(tickData);
                    tickData = [];
                    zticks(tickData);
                end
                
                tightInset = get(gca, 'TightInset');
                set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                           (tightInset(2) + 0.00625), ...
                                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
                pause(0.5);
                hold off;
                
                % Save Figure
                print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
                
                if figSave
                    savefig(gcf, [userpath, '/Output/Figures/', figName, '.fig']);
                end
                
            else
                planeNo = planeNo + 1;
            end
    
    end
    
end
