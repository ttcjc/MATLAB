%% Lagrangian Volume Field Plotter v1.0
% ---
% Plots Previously Processed Contaminant Deposition Maps
% Usage:
% ---


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function fig = volumeFieldPlots(x, y, z, fieldData, ...
                                fig, figName, geometry, isoValue, ...
                                figTitle, xLimsPlot, yLimsPlot, zLimsPlot)

    fieldData = smooth3(fieldData, 'gaussian');
    
    % Figure Setup
    fig = fig + 1;
    set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
    set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
             'fontSize', 18, 'layer', 'top');
    hold on;
    
    % Plot
    part = fieldnames(geometry);
    for j = 1:height(part)
        patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                    'edgeColor', [0.5, 0.5, 0.5], ...
                                    'lineStyle', 'none');
    end
    
    surface = isosurface(x, y, z, fieldData, isoValue);
    patch(surface, 'faceColor', [0.21176, 0.06667, 0.38824], 'edgeColor', 'none');
    
    % Figure Formatting
    title(figTitle);
    lightangle(0, 45);
    lighting gouraud;
    box on;
    view([30, 30]);
    xlim([xLimsPlot(1), xLimsPlot(2)]);
    ylim([yLimsPlot(1), yLimsPlot(2)]);
    zlim([zLimsPlot(1), zLimsPlot(2)]);
    xTickData = xLimsPlot(1):((xLimsPlot(2) - xLimsPlot(1)) / 5):xLimsPlot(2);
    xticks(xTickData(2:(end-1)));
    yTickData = yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2);
    yticks(yTickData(2:(end-1)));
    zTickData = zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2);
    zticks(zTickData(2:(end-1)));
    xtickformat('%+.3f');
    ytickformat('%+.3f');
    ztickformat('%+.3f');
    xT = xlabel('x_{\it{l}}');
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