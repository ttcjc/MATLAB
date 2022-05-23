%% Lagrangian Contaminant Map Plotter v2.0
% ----
% Plots Previously Processed Contaminant Maps
% ----
% Usage: fig = vectorPlots();


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Main Function

function fig = contaminantPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                                basePerim, positionData, contaminantData, fig, figName, cMap, geometry, ...
                                xDims, CoM, figTitle, figSubtitle, cLims, normalise)
    
    cellSize = 0.5e-3; % [m or l]
    
    % Generate Gridded Data
    cellSizeX = cellSize;
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));
    
    [x, y, z] = meshgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                         yLimsData(1):cellSizeY:yLimsData(2), ...
                         zLimsData(1):cellSizeZ:zLimsData(2));
    
    interp = scatteredInterpolant(positionData(:,2), positionData(:,3), contaminantData(:,1), ...
                                  'linear', 'none');
    
    contamination = zeros(size(x));
    contamination(:,2,:) = interp(y(:,2,:), z(:,2,:));
    
    if ~isempty(basePerim)
        [indexIn, indexOn] = inpolygon(y, z, basePerim(:,2), basePerim(:,3));
        indexBase = double(or(indexIn, indexOn));
        indexBase(indexBase == 0) = nan;

        contamination = contamination .* indexBase;
    end

    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    colormap(cMap);
    hold on;
    
    % Figure Plotting
    parts = fieldnames(geometry);
    
    for i = 1:height(parts)
        patch('faces', geometry.(parts{i}).faces, ...
              'vertices', geometry.(parts{i}).vertices, ...
              'faceColor', ([128, 128, 128] / 255), ...
              'edgeColor', ([128, 128, 128] / 255), ...
              'lineStyle', 'none');
    end
    
    surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(contamination(:,2,:)), ...
         'lineStyle', 'none', 'faceLighting', 'none');

    if xLimsData < xDims(1) || xLimsData > xDims(2)
        
        for i = 1:height(parts)
            geometry.(parts{i}).boundaries.YZ(:,1) = xLimsData;

            plot3(geometry.(parts{i}).boundaries.YZ(:,1), ...
                  geometry.(parts{i}).boundaries.YZ(:,2), ...
                  geometry.(parts{i}).boundaries.YZ(:,3), ...
                  'color', 'w', 'lineStyle', '-', 'lineWidth', 1.5);
        end
        
    end
    
    scatter3(CoM(1), CoM(2), CoM(3), 125, 'w', 'filled')

    % Figure Formatting
    title(figTitle, 'color', ([254, 254, 254] / 255));
    subtitle(figSubtitle);
    lightangle(90, 45);
    axis on;
    box on;
    caxis(cLims);
    view([90, 0]);
    xlim([xLimsPlot(1), xLimsPlot(2)]);
    ylim([yLimsPlot(1), yLimsPlot(2)]);
    zlim([zLimsPlot(1), zLimsPlot(2)]);
    tickData = [];
    xticks(tickData);
    tickData = round((yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2)), 2);
    yticks(tickData(2:end-1));
    tickData = round((zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2)), 2);
    zticks(tickData(2:end-1));
    xtickformat('%+.2f');
    ytickformat('%+.2f');
    ztickformat('%+.2f');

    if normalise
        xT = xlabel([]);
        yT = ylabel('y_{\it{l}}');
        zT = zlabel('z_{\it{l}}');
    else
        xT = xlabel([]);
        yT = ylabel('y_{\it{m}}');
        zT = zlabel('z_{\it{m}}');
    end

    xT.FontName = 'LM Roman 12';
    yT.FontName = 'LM Roman 12';
    zT.FontName = 'LM Roman 12';
    hold off;

    pause(2);
    exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
    
end