%% Volume Field Plotter v2.0
% ----
% Plots Previously Processed Volume Fields
% ----
% Usage: fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, ...
%                               fieldDataA, fieldDataB, fig, figName, geometry, isoValue, ...
%                               cMap, fieldColour, figTitle, figSubtitle, ...
%                               xLimsPlot, yLimsPlot, zLimsPlot);


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Main Function

function fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, ...
                                fieldDataA, fieldDataB, fig, figName, geometry, isoValue, ...
                                cMap, fieldColour, figTitle, figSubtitle, ...
                                xLimsPlot, yLimsPlot, zLimsPlot)
                            
    % Generate Refined Grid
    cellSize = 4e-3;

    cellSizeX = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));
    
    [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                       yLimsData(1):cellSizeY:yLimsData(2), ...
                       zLimsData(1):cellSizeZ:zLimsData(2));

    % Smooth Data
    if POD
        fieldDataA = interpn(xInit, yInit, zInit, fieldDataA, x, y, z);
        fieldDataB = interpn(xInit, yInit, zInit, fieldDataB, x, y, z);
        
        fieldDataA = smooth3(fieldDataA, 'gaussian');
        fieldDataB = smooth3(fieldDataB, 'gaussian');
    else
        fieldDataA = interpn(xInit, yInit, zInit, fieldDataA, x, y, z);
        
        fieldDataA = smooth3(fieldDataA, 'gaussian');
    end
    
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
        iso = isosurface(x, y, z, fieldDataA, isoValue);
        iso = patch(iso, 'faceColor', cMap(1,:), 'edgeColor', 'none');
        isonormals(fieldDataA, iso);
        
        iso = isosurface(x, y, z, fieldDataB, isoValue);
        iso = patch(iso, 'faceColor', cMap(end,:), 'edgeColor', 'none');
        isonormals(fieldDataB, iso);
    else
        iso = isosurface(x, y, z, fieldDataA, isoValue);
        iso = patch(iso, 'faceColor', fieldColour, 'edgeColor', 'none');
        isonormals(fieldDataA, iso);
    end
    
    % Figure Formatting
    title(figTitle, 'color', ([254, 254, 254] / 255));
    subtitle(figSubtitle);
    lightangle(0, 45);
    axis on;
    box on;
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
    set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
    hold off;
    
    pause(2);
    exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);