clc;

close all;

test = unique(cat(1, LagDataPlane.X_P18637.positionCartesian{:}), 'rows');

disp('Select an Impact Point (Right-Click):')

valid = false;
while ~valid
    set(figure(1), 'name', 'Impacts', 'color', [1, 1, 1], 'units', 'pixels');
    set(gca, 'dataAspectRatio', [1, 1, 1], 'lineWidth', 4, 'fontName', 'LM Mono 12', ...
             'fontSize', 18, 'layer', 'top');
    hold on;

    cMap = viridis(3); cMap = cMap(1,:);

    scatter(test(:,2), test(:,3), 5, cMap, 'filled');

    xlim([-ceil(max(abs(test(:,2)))), ceil(max(abs(test(:,2))))]);
    ylim([0, ceil(max(test(:,3)))]);
    
    try
        [xi,yi] = getpts;
        
        if numel(xi) > 1
            disp('    Warning: Multiple Impact Positions Selected')
        else
            close(1);
            valid = true;
        end
        
    catch
        disp('    Warning: A Valid Impact Position Must Be Selected')
    end
    
end