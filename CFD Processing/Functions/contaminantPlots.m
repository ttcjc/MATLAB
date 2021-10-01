%% Lagrangian Contaminant Deposition Map Plotter v1.0
% ---
% Plots Previously Processed Contaminant Deposition Maps
% Usage: fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
%                               positionData, contaminantData, ...
%                               fig, figName, cMap, geometry, ...
%                               CoM, modelOutline, figTitle, cLims, ...
%                               xLimsPlot, yLimsPlot, zLimsPlot)
%        '*LimsData' -> Contour Plot Limits
%        'positionData' -> Cartesian Positions
%        'contaminantData' -> Deposition Variable @ 'positionData'
%        'fig' -> Figure Number
%        'figName' -> Figure Name
%        'cMap' -> Colour Map
%        'geometry' -> STL to Include in Plot
%        'CoM' -> Cartesian Position of Centre of Mass
%        'modelOutline' -> Cartesian Positions of Model Vertices
%        'figTitle' -> Figure Title
%        'cLims' -> Colour Bounds
%        '*LimsPlot' -> Figure Limits
% ---


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = contaminantPlots(xLimsData, yLimsData, zLimsData, ...
                                positionData, contaminantData, ...
                                fig, figName, cMap, geometry, ...
                                CoM, modelOutline, figTitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot)
    
    % Generate Gridded Data
    cellSize = 0.001; % l
    
    cellSizeX = cellSize;
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round((yLimsData(2) - yLimsData(1)) / cellSize);
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round((zLimsData(2) - zLimsData(1)) / cellSize);
    
    [x, y, z] = meshgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                         yLimsData(1):cellSizeY:yLimsData(2), ...
                         zLimsData(1):cellSizeZ:zLimsData(2));

    interp = scatteredInterpolant(positionData(:,2), ...
                                  positionData(:,3), ...
                                  contaminantData(:,1), ...
                                  'linear', 'linear');

    deposition = zeros(size(x));
    deposition(:,2,:) = interp(y(:,2,:), z(:,2,:));

    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
             'fontSize', 18, 'layer', 'top');
    colormap(cMap);
    hold on;

    % Plot
    part = fieldnames(geometry);
    for j = 1:height(part)
        patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                    'edgeColor', [0.5, 0.5, 0.5], ...
                                    'lineStyle', 'none');
    end

    surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(deposition(:,2,:)), ...
         'lineStyle', 'none', 'faceLighting', 'none');

    if ~isempty(CoM)
        scatter3(CoM(1,1), CoM(1,2), CoM(1,3), 125, 'w', '*', 'lineWidth', 1.25);
    end
     
    if ~isempty(modelOutline)
        plot3(modelOutline(:,1), modelOutline(:,2), modelOutline(:,3), ...
              'color', 'w', 'lineStyle', ':', 'lineWidth', 1.25);
    end

    % Figure Formatting
    title(figTitle);
    lightangle(90, 45);
    box on;
    caxis(cLims);
    view([90, 0]);
    xlim([xLimsPlot(1), xLimsPlot(2)]);
    ylim([yLimsPlot(1), yLimsPlot(2)]);
    zlim([zLimsPlot(1), zLimsPlot(2)]);
    xTickData = [];
    xticks(xTickData);
    yTickData = yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2);
    yticks(yTickData(2:(end-1)));
    zTickData = zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2);
    zticks(zTickData(2:(end-1)));
    xtickformat('%+.3f');
    ytickformat('%+.3f');
    ztickformat('%+.3f');
    xT = xlabel([]);
    yT = ylabel('y_{\it{l}}');
    zT = zlabel('z_{\it{l}}');
    xT.FontName = 'LM Roman 12';
    yT.FontName = 'LM Roman 12';
    zT.FontName = 'LM Roman 12';
    hold off;

    % Save Figure
    pause(2);
    exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
%     savefig(fig, ['~/MATLAB/Output/Figures/', figName]);
    
end