%% Volume Field Plotter v2.0
% ----
% Plots Previously Processed Volume Fields
% ----
% Usage: fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
%                               fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
%                               xLimsPlot, yLimsPlot, zLimsPlot);
%        '*LimsData'   -> Contour Plot Limits
%        '*Init'       -> Initial 3D Arrays of Cartesian Positions
%        'POD'         -> POD Mode Presentation [True/False]
%        'fieldData'   -> 3D Array of Field Data @ '*Init' Points
%        'fig'         -> Figure Number
%        'figName'     -> Figure Name
%        'geometry'    -> STL(s) to Include in Plot
%        'isoValue'    -> Field Value Used for Isosurface
%        'cMap'        -> Colour Map
%        'figTitle'    -> Leave Blank ('-') for Formatting Purposes
%        'figSubtitle' -> Figure Title
%        '*LimsPlot'   -> 3D Axes Limits


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Cleaned-up POD Functionality


%% Main Function

function fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                                fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                                xLimsPlot, yLimsPlot, zLimsPlot)
    
    % Generate Refined Grid
    cellSize = 4e-3;

    cellSizeX = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize));
    cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize));
    cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize));
    
    [x, y, z] = ndgrid(xLimsData(1):cellSizeX:xLimsData(2), ...
                       yLimsData(1):cellSizeY:yLimsData(2), ...
                       zLimsData(1):cellSizeZ:zLimsData(2));
    
    % Convert From 'ndgrid' to 'meshgrid' Format
    xInit = permute(xInit, [2,1,3]);
    yInit = permute(yInit, [2,1,3]);
    zInit = permute(zInit, [2,1,3]);
    x = permute(x, [2,1,3]);
    y = permute(y, [2,1,3]);
    z = permute(z, [2,1,3]);

    % Smooth Data
    if POD
        
        if iscell(fieldData) && height(fieldData) == 2
            
            for i = 1:height(fieldData)
                fieldData{i} = permute(fieldData{i}, [2,1,3]);
                fieldData{i} = interp3(xInit, yInit, zInit, fieldData{i}, x, y, z);
                fieldData{i} = smooth3(fieldData{i}, 'gaussian');
            end
            
        else
            error('Presentation of POD Modes Necessitates ''fieldData'' Be in the Form of a {2,1} Cell Array');
        end
        
    else
        
        if ~iscell(fieldData)
            fieldData = permute(fieldData, [2,1,3]);
            fieldData = interp3(xInit, yInit, zInit, fieldData, x, y, z);
            fieldData = smooth3(fieldData, 'gaussian');
        else
            error('''fieldData'' Must Be a 3D Array of Type ''double''');
        end
    end
    
    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'color', [1, 1, 1], 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    lighting gouraud;
    hold on;
    
    % Figure Plotting
    if ~isempty(geometry)
        parts = fieldnames(geometry);

        for i = 1:height(parts)
            patch('faces', geometry.(parts{i}).faces, ...
                  'vertices', geometry.(parts{i}).vertices, ...
                  'faceColor', [0.5, 0.5, 0.5], ...
                  'edgeColor', [0.5, 0.5, 0.5], ...
                  'lineStyle', 'none');
        end

    end

    if POD
        iso = isosurface(x, y, z, fieldData{1}, isoValue);
        iso = patch(iso, 'faceColor', cMap(1,:), 'edgeColor', 'none');
        iso.VertexNormals = isonormals(x, y, z, fieldData{1}, iso);
        
        iso = isosurface(x, y, z, fieldData{2}, isoValue);
        iso = patch(iso, 'faceColor', cMap(end,:), 'edgeColor', 'none');
        iso.VertexNormals = isonormals(x, y, z, fieldData{2}, iso);
    else
        iso = isosurface(x, y, z, fieldData, isoValue);
        iso = patch(iso, 'faceColor', cMap(1,:), 'edgeColor', 'none');
        iso.VertexNormals = isonormals(x, y, z, fieldData, iso);
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
    exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
    savefig(gcf, ['~/MATLAB/Output/Figures/', figName, '.fig']);

end