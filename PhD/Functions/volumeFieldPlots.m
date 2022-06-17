%% Volume Field Plotter v2.0
% ----
% Plots Previously Processed Volume Fields
% ----
% Usage: fig = volumeFieldPlots()


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Main Function

function fig = volumeFieldPlots(xLims, yLims, zLims, volumeData, fig, figName, geometry, POD, xInit, yInit, zInit, ...
                                isoValue, figTitle, figSubtitle)

    % Generate Refined Grid
    cellSize = 4e-3;

    cellSizeX = (yLims(2) - yLims(1)) / round(((yLims(2) - yLims(1)) / cellSize));
    cellSizeY = (yLims(2) - yLims(1)) / round(((yLims(2) - yLims(1)) / cellSize));
    cellSizeZ = (zLims(2) - zLims(1)) / round(((zLims(2) - zLims(1)) / cellSize));
    
    [x, y, z] = ndgrid(xLims(1):cellSizeX:xLims(2), ...
                       yLims(1):cellSizeY:yLims(2), ...
                       zLims(1):cellSizeZ:zLims(2));

    volumeData = interp3(xInit, yInit, zInit, volumeData, x, y, z);

    % Smooth Data
    volumeData = smooth3(volumeData, 'gaussian');
    
    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    lighting gouraud;
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

    if POD
        iso = isosurface(x, y, z, volumeData, isoValue);
        iso = patch(iso, 'faceColor', max(cMap), 'edgeColor', 'none');
        isonormals(volumData, iso);
        
        iso = isosurface(x, y, z, volumeData, -isoValue);
        iso = patch(iso, 'faceColor', min(cMap), 'edgeColor', 'none');
        isonormals(volumData, iso);
    else
        iso = isosurface(x, y, z, volumeData, isoValue);
        iso = patch(iso, 'faceColor', fieldColour, 'edgeColor', 'none');
        isonormals(volumData, iso);
    end
    
    % Figure Formatting
    title(figTitle, 'color', ([254, 254, 254] / 255));
    subtitle(figSubtitle);
    lightangle(0, 45);
    axis on;
    box on;
    view([30, 30]);
    xlim([xLims(1), xLims(2)]);
    ylim([yLims(1), yLims(2)]);
    zlim([zLims(1), zLims(2)]);
    tickData = [];
    xticks(tickData);
    tickData = [];
    yticks(tickData);
    tickData = [];
    zticks(tickData);
    set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
    hold off;
    
    pause(2);
    exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);