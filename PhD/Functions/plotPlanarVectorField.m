%% Planar Vector Field Plotter v2.4
% ----
% Plots Previously Processed Planar Vector Fields
% ----
% Usage: [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
%                                               xLimsData, yLimsData, zLimsData, nComponents, component, ...
%                                               mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
%                                               streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
%                                               zLimsPlot, normDims, figSave);
%
%        'orientation'  -> Plane Orientation ['YZ', 'XZ', 'XY']
%        'positionData' -> Cartesian Positions of Data Points
%        'vectorData'   -> Field Values @ 'positionData' Points
%        'spatialRes'   -> Target Grid Spacing [Dimensions of 'positionData']
%        '*LimsData'    -> Contour Plot Limits [Dimensions of 'positionData']
%        'nComponents'  -> Number of Vector Components Used to Calculate Magnitude
%        'component'    -> Specific Component Used When 'nComponents' Is Set to 1
%        'mapPerim'     -> Map Perimeter Used When Plotting a Plane of Arbitrary Shape [Dimensions of 'positionData']
%        'nPlanes'      -> Number of Planes in a Multi-Plane Figure
%        'planeNo'      -> Current Plane Number
%        'fig'          -> Figure Number
%        'figName'      -> Figure Name
%        'cMap'         -> Colour Map
%        'geometry'     -> STL(s) to Include in Plot
%        'streamlines'  -> Include Streamlines [True/False]
%        'figTitle'     -> Figure Title
%        'cLims'        -> Colour Map Limits
%        '*LimsPlot'    -> 3D Axes Limits [Dimensions of 'positionData']
%        'normDims'     -> Normalise Dimensions [True/False]
%        'figSave'      -> Save .fig File [True/False]


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Shifted to Using 'griddedInterpolant' and Added Support for Multi-Plane Figures)
% v2.1 - Rename and Minor Formatting Updates
% v2.2 - Added Spatial Resolution as an Input Variable
% v2.3 - Update To Ensure Consistent Figure Sizes
% v2.4 - Improved Consistency of Grid Spacing


%% Main Function

function [fig, planeNo] = plotPlanarVectorField(orientation, positionData, vectorData, spatialRes, ...
                                                xLimsData, yLimsData, zLimsData, nComponents, component, ...
                                                mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, ...
                                                streamlines, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                zLimsPlot, normDims, figSave)
    
    % Format Data
    switch orientation
        
        case 'YZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,2))), ...
                         height(unique(positionData(:,3)))];

            y = reshape(positionData(:,2), gridShape);
            z = reshape(positionData(:,3), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(y, z, u, 'linear', 'none');
            vInterp = griddedInterpolant(y, z, v, 'linear', 'none');
            wInterp = griddedInterpolant(y, z, w, 'linear', 'none');

            % Generate 3D Surface
            nPy = (diff(yLimsData) / spatialRes) + 1;
            nPz = (diff(zLimsData) / spatialRes) + 1;
            
            [x, y, z] = ndgrid((xLimsData - spatialRes):spatialRes:(xLimsData + spatialRes), ...
                               linspace(yLimsData(1), yLimsData(2), nPy), ...
                               linspace(zLimsData(1), zLimsData(2), nPz));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(2,:,:) = uInterp(y(2,:,:), z(2,:,:));

            v = zeros(size(x));
            v(2,:,:) = vInterp(y(2,:,:), z(2,:,:));

            w = zeros(size(x));
            w(2,:,:) = wInterp(y(2,:,:), z(2,:,:));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(2,:,:);
            elseif nComponents == 2
                vector = sqrt(v(2,:,:).^2 + w(2,:,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(2,:,:).^2 + v(2,:,:).^2 + w(2,:,:).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(y, z, mapPerim(:,2), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
            end

        case 'XZ'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,3)))];

            x = reshape(positionData(:,1), gridShape);
            z = reshape(positionData(:,3), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(x, z, u, 'linear', 'none');
            vInterp = griddedInterpolant(x, z, v, 'linear', 'none');
            wInterp = griddedInterpolant(x, z, w, 'linear', 'none');

            % Generate 3D Surface
            nPx = (diff(xLimsData) / spatialRes) + 1;
            nPz = (diff(zLimsData) / spatialRes) + 1;

            [x, y, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                               (yLimsData - spatialRes):spatialRes:(yLimsData + spatialRes), ...
                               linspace(zLimsData(1), zLimsData(2), nPz));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(:,2,:) = uInterp(x(:,2,:), z(:,2,:));

            v = zeros(size(x));
            v(:,2,:) = vInterp(x(:,2,:), z(:,2,:));

            w = zeros(size(x));
            w(:,2,:) = wInterp(x(:,2,:), z(:,2,:));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,2,:);
            elseif nComponents == 2
                vector = sqrt(u(:,2,:).^2 + w(:,2,:).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,2,:).^2 + v(:,2,:).^2 + w(:,2,:).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, z, mapPerim(:,1), mapPerim(:,3));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
            end
            
        case 'XY'
            % Reshape Data for Improved Interpolation Performance
            gridShape = [height(unique(positionData(:,1))), ...
                         height(unique(positionData(:,2)))];

            x = reshape(positionData(:,1), gridShape);
            y = reshape(positionData(:,2), gridShape);

            u = reshape(vectorData(:,1), gridShape);
            v = reshape(vectorData(:,2), gridShape);
            w = reshape(vectorData(:,3), gridShape);

            % Perform Interpolation
            uInterp = griddedInterpolant(x, y, u, 'linear', 'none');
            vInterp = griddedInterpolant(x, y, v, 'linear', 'none');
            wInterp = griddedInterpolant(x, y, w, 'linear', 'none');

            % Generate 3D Surface
            nPx = (diff(xLimsData) / spatialRes) + 1;
            nPy = (diff(yLimsData) / spatialRes) + 1;
            
            [x, y, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                               linspace(yLimsData(1), yLimsData(2), nPy), ...
                               (zLimsData - spatialRes):spatialRes:(zLimsData + spatialRes));
            
            % Map Data on to 3D Surface
            u = zeros(size(x));
            u(:,:,2) = uInterp(x(:,:,2), y(:,:,2));

            v = zeros(size(x));
            v(:,:,2) = vInterp(x(:,:,2), y(:,:,2));

            w = zeros(size(x));
            w(:,:,2) = wInterp(x(:,:,2), y(:,:,2));

            % Evaluate Vector Components for Contour
            if nComponents == 1
                vector = eval(component);
                vector = vector(:,:,2);
            elseif nComponents == 2
                vector = sqrt(u(:,:,2).^2 + v(:,:,2).^2);
            elseif nComponents == 3
                vector = sqrt(u(:,:,2).^2 + v(:,:,2).^2 + w(:,:,2).^2);
            end

            % Remove Data Outside Desired Map Perimeter
            if ~isempty(mapPerim)
                [indexIn, indexOn] = inpolygon(x, y, mapPerim(:,1), mapPerim(:,2));
                indexMap = double(or(indexIn, indexOn));

                vector = vector .* indexMap;
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
                 'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
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
            % Plot Vector Field
            surf(squeeze(x(2,:,:)), squeeze(y(2,:,:)), squeeze(z(2,:,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);

            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                v = permute(v, [2,1,3]); v(isnan(v)) = 0;
                w = permute(w, [2,1,3]); w(isnan(w)) = 0;
                
                streams = streamslice(x, y, z, zeros(size(x)), v, w, xLimsData, [], [], ...
                                      1.25, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
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
                    tickData = (yLimsPlot(1):(diff(yLimsPlot) / 5):yLimsPlot(2));
                    yticks(tickData(2:5));
                    tickData = (zLimsPlot(1):(diff(zLimsPlot) / 5):zLimsPlot(2));
                    zticks(tickData(2:5));
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normDims
                        ylabel({'{$y_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
                        zlabel({'{-----}'; '{$z_{\ell}$}'}, 'interpreter', 'latex');
                    else
                        ylabel({'{$y$ ($m$)}'; '{-----}'}, 'interpreter', 'latex');
                        zlabel({'{-----}'; '{$z$ ($m$)}'}, 'interpreter', 'latex');
                    end

                else
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
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
            % Plot Vector Field
            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);

            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]); u(isnan(u)) = 0;
                w = permute(w, [2,1,3]); w(isnan(w)) = 0;
                
                streams = streamslice(x, y, z, u, zeros(size(x)), w, [], yLimsData, [], ...
                                      1.25, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([0, 0]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):(diff(xLimsPlot) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:5));
                    tickData = [];
                    yticks(tickData);
                    tickData = round((zLimsPlot(1):(diff(zLimsPlot) / 5):zLimsPlot(2)), 2);
                    zticks(tickData(2:5));
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normDims
                        xlabel({'{$x_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
                        zlabel({'{-----}'; '{$z_{\ell}$}'}, 'interpreter', 'latex');
                    else
                        xlabel({'{$x$ ($m$)}'; '{-----}'}, 'interpreter', 'latex');
                        zlabel({'{-----}'; '{$z$ ($m$)}'}, 'interpreter', 'latex');
                    end

                else
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
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
            % Plot Vector Field
            surf(squeeze(x(:,:,2)), squeeze(y(:,:,2)), squeeze(z(:,:,2)), squeeze(vector), ...
                 'lineStyle', 'none', 'faceLighting', 'none', 'faceAlpha', 0.95);
            
            % Plot Streamlines
            if streamlines
                % Convert From 'ndgrid' to 'meshgrid' Format
                x = permute(x, [2,1,3]);
                y = permute(y, [2,1,3]);
                z = permute(z, [2,1,3]);
                u = permute(u, [2,1,3]); u(isnan(u)) = 0;
                v = permute(v, [2,1,3]); v(isnan(v)) = 0;
                
                streams = streamslice(x, y, z, u, v, zeros(size(x)), [], [], zLimsData, ...
                                      1.25, 'arrows', 'linear');
                set(streams, 'color', 'k', 'lineStyle', '-', 'lineWidth', 2);
            end
            
            % Format Figure
            if planeNo == nPlanes

                if nPlanes == 1
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
                    lightangle(0, 45);
                    axis on;
                    box on;
                    grid off;
                    caxis(cLims);
                    view([0, 90]);
                    xlim([xLimsPlot(1), xLimsPlot(2)]);
                    ylim([yLimsPlot(1), yLimsPlot(2)]);
                    zlim([zLimsPlot(1), zLimsPlot(2)]);
                    tickData = round((xLimsPlot(1):(diff(xLimsPlot) / 5):xLimsPlot(2)), 2);
                    xticks(tickData(2:5));
                    tickData = round((yLimsPlot(1):(diff(yLimsPlot) / 5):yLimsPlot(2)), 2);
                    yticks(tickData(2:5));
                    tickData = [];
                    zticks(tickData);
                    xtickformat('%+.2g');
                    ytickformat('%+.2g');
                    ztickformat('%+.2g');
                    
                    if normDims
                        xlabel({'{$x_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
                        ylabel({'{-----}'; '{$y_{\ell}$}'}, 'interpreter', 'latex');
                    else
                        xlabel({'{$x$ ($m$)}'; '{-----}'}, 'interpreter', 'latex');
                        ylabel({'{-----}'; '{$y$ ($m$)}'}, 'interpreter', 'latex');
                    end

                else
                    title('{-----}', 'interpreter', 'latex');
                    subtitle(figTitle);
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