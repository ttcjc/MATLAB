%% Numerical Pressure Processing v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('==================================');
disp ('Numerical Pressure Processing v1.0');
disp ('==================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Initialisation

[caseFolder, pressureData, xDims, yDims, zDims, geometry] = initialisePVdata('p');

disp(' ');
disp(' ');

% Temporary Implementation Control
if ~contains(caseFolder, 'Windsor')
    error('Case Type Not Yet Supported');
end


%% Process Pressure Data

disp('PROCESS PRESSURE DATA');
disp('---------------------');
disp(' ');
disp('***********');
disp('  Running  ');

tic;
              
% Set and Normalise Dimensions
xPre = max(width(extractAfter(num2str(xDims(1), 8), '.')), width(extractAfter(num2str(xDims(2), 8), '.')));
yPre = max(width(extractAfter(num2str(yDims(1), 8), '.')), width(extractAfter(num2str(yDims(2), 8), '.')));
zPre = max(width(extractAfter(num2str(zDims(1), 8), '.')), width(extractAfter(num2str(zDims(2), 8), '.')));

if contains(caseFolder, 'Windsor')
    
    if contains(caseFolder, 'Balance')
        Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
        At = (2 * (0.96 + (4.704 * tan(atan(0.02613 / 9.408)))) * 1.32);
    elseif contains(caseFolder, 'Upstream')
        Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
        At = (2 * (0.96 + (3.379 * tan(atan(0.02613 / 9.408)))) * 1.32);
        xDims = xDims + 1.325;
    end
    
    xDims = round(xDims / 1.044, xPre);
    yDims = round(yDims / 1.044, yPre);
    zDims = round(zDims / 1.044, zPre);    
else
    % Add Support for Future Geometries
end

% Identify Model Boundaries
part = fieldnames(geometry);
for i = 1:height(part)

    if round(max(geometry.(part{i,1}).vertices(:,1)), xPre) == xDims(2)
        break
    end

    if i == height(part)
        disp(' ');
        error('Mismatch Between ''xDims'' and Geometry Bounding Box')
    end

end

geoPoints = unique(geometry.(part{i,1}).vertices, 'stable', 'rows');
index = find(round(geoPoints(:,1), xPre) == xDims(2));

xDimsBase = xDims(2);
yDimsBase = round([min(geoPoints(index,2)); max(geoPoints(index,2))], yPre);
zDimsBase = round([min(geoPoints(index,3)); max(geoPoints(index,3))], zPre);

disp(' ');

% Format Data
for i = 1:height(pressureData.files)
    plane = pressureData.files{i,1}((max(strfind(pressureData.files{i,1}, '_')) + 1):(end - 4));
    disp(['    Formatting ', plane, ' Data']);
    
    import = readmatrix([caseFolder, '/', pressureData.files{i,1}]);
    
    if contains(caseFolder, 'Windsor')
        U = 40; % m/s
        rho = 1.269; % kg/m^3
        pRef = 0 * rho; % Pa
        
        if contains(caseFolder, 'Upstream')
            import(:,1) = import(:,1) + 1.325;
        end
        
        switch plane

            case 'Base'
                import = import((round(import(:,1), xPre) == 0.48325),:);
                
                xDimsBase = xDimsBase + 0.002;
                yDimsBase(1) = yDimsBase(1) + 0.002;
                yDimsBase(2) = yDimsBase(2) - 0.002;
                zDimsBase(1) = zDimsBase(1) + 0.002;
                zDimsBase(2) = zDimsBase(2) - 0.002;
        
                xLims = round([-0.61075; 0.53325] / 1.044, xPre);
                yLims = round([-0.2445; 0.2445] / 1.044, yPre);
                zLims = round([0; 0.389] / 1.044, zPre);
                
            case 'Centreline'
                import = import((round(import(:,3), zPre) > 0.125),:);
                import = sortrows(import, 1);
                
                xLims = [(xDims(1) - 0.02), (xDims(2) + 0.02)];

        end
    
    else
        % Add Support for Future Geometries
    end
    
    pressureData.(plane).position(:,1) = round(import(:,1) / 1.044, xPre);
    pressureData.(plane).position(:,2) = round(import(:,2) / 1.044, yPre);
    pressureData.(plane).position(:,3) = round(import(:,3) / 1.044, zPre);
    pressureData.(plane).p = import(:,4) * rho;
    pressureData.(plane).Cp = (pressureData.(plane).p - pRef) / (0.5 * rho * U^2);
    
    % Remove Duplicate Entries
    [pressureData.(plane).position, index, ~] = unique(pressureData.(plane).position, ...
                                                       'rows', 'stable');
    pressureData.(plane).p = pressureData.(plane).p(index,:);
    pressureData.(plane).Cp = pressureData.(plane).Cp(index,:);   
    
    % Calculate Centre of Pressure
    switch plane
        
        case 'Base'
            pressureData.(plane).CoP = zeros(1,3);
            pressureData.(plane).CoP(1,1) = xDimsBase;
            pressureData.(plane).CoP(1,2) = sum(pressureData.(plane).position(:,2) .* pressureData.(plane).p) / ...
                                            sum(pressureData.(plane).p);
            pressureData.(plane).CoP(1,3) = sum(pressureData.(plane).position(:,3) .* pressureData.(plane).p) / ...
                                            sum(pressureData.(plane).p);
                                        
    end
    
    % Blockage Correction
    E = Am / At;
    pressureData.(plane).CpCorr = (pressureData.(plane).Cp + (2 * E)) / (1 + (2 * E));
end

disp(' ')

for i = 1:height(pressureData.files)
    plane = pressureData.files{i,1}((max(strfind(pressureData.files{i,1}, '_')) + 1):(end - 4));
    disp(['    Presenting ', plane, ' Data']);
    
    switch plane
        
        case 'Base'
            cellSize = 0.001; % l
            cellSizeX = cellSize;
            cellSizeY = (2 * yDimsBase(2)) / round((2 * yDimsBase(2)) / cellSize);
            cellSizeZ = (zDimsBase(2) - zDimsBase(1)) / round((zDimsBase(2) - zDimsBase(1)) / cellSize);
            
            [x, y, z] = meshgrid((xDimsBase - cellSizeX):cellSizeX:(xDimsBase + cellSizeX), ...
                         yDimsBase(1):cellSizeY:yDimsBase(2), ...
                         zDimsBase(1):cellSizeZ:zDimsBase(2));

            interp = scatteredInterpolant(pressureData.(plane).position(:,2), ...
                                          pressureData.(plane).position(:,3), ...
                                          pressureData.(plane).CpCorr(:,1), ...
                                          'linear', 'none');
                                      
            Cp = zeros(size(x));
            Cp(:,2,:) = interp(y(:,2,:), z(:,2,:));
            
            % Figure Setup
            fig = fig + 1;
            figName = ['Cp_', plane, '_CFD'];
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [1, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 18, 'layer', 'top');
            colormap(viridis(24));
            hold on;

            % Plot
            part = fieldnames(geometry);
            for j = 1:height(part)
                patch(geometry.(part{j,1}), 'faceColor', [0.5, 0.5, 0.5], ...
                                            'edgeColor', [0.5, 0.5, 0.5], ...
                                            'lineStyle', 'none');
            end

            surf(squeeze(x(:,2,:)), squeeze(y(:,2,:)), squeeze(z(:,2,:)), squeeze(Cp(:,2,:)), ...
                 'lineStyle', 'none', 'faceLighting', 'none');
             
            scatter3(pressureData.(plane).CoP(1,1), pressureData.(plane).CoP(1,2), pressureData.(plane).CoP(1,3), 125, 'w', '*', 'lineWidth', 1.25);

            % Figure Formatting
            title(' ');
            light;
            box on;
            caxis([-0.24, -0.06]);
            % caxis([-0.226860044655411, -0.075288561534909]); % Square-Back CFD
            % caxis([-0.235435242242635, -0.057290950577478]); % Side Tapers CFD
            % caxis([-0.234098930226398, -0.127384500508930]); % Square-Back Exp
            % caxis([-0.210908693597759, -0.098457129729324]); % Side Tapers Exp
%             caxis([-0.235435242242635, -0.057290950577478]); % Comparison
            view([90, 0]);
            xlim([xLims(1), xLims(2)]);
            ylim([yLims(1), yLims(2)]);
            zlim([zLims(1), zLims(2)]);
            tickData = [];
            xticks(tickData);
            tickData = yLims(1):((yLims(2) - yLims(1)) / 5):yLims(2);
            yticks(tickData(2:(end-1)));
            tickData = zLims(1):((zLims(2) - zLims(1)) / 5):zLims(2);
            zticks(tickData(2:(end-1)));
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

            pause(2);
            exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
            
        case 'Centreline'
            % Figure Setup
            fig = fig + 1;
            figName = ['Cp_', plane, '_CFD'];
            set(figure(fig), 'outerPosition', [25, 25, 850, 850], 'name', figName);
            set(gca, 'dataAspectRatio', [0.25, 1, 1], 'fontName', 'LM Mono 12', ...
                     'fontSize', 18, 'layer', 'top');
            hold on;
            
            % Plot
            plot(pressureData.(plane).position(:,1), pressureData.(plane).CpCorr, ...
                'color', [0.21176, 0.06667, 0.38824]);
            
            % Figure Formatting
            title(' ');
            box on;
            grid on;
            xlim([xLims(1), xLims(2)]);
            ylim([-1, 1]);
            tickData = xLims(1):((xLims(2) - xLims(1)) / 5):xLims(2);
            xticks(tickData(2:(end-1)));
            tickData = -1:0.4:1;
            yticks(tickData(2:(end-1)));
            xtickformat('%+.3f');
            ytickformat('%+.3f');
            xT = xlabel('x_{\it{l}}');
            yT = ylabel('C_p');
            xT.FontName = 'LM Roman 12';
            yT.FontName = 'LM Roman 12';
            hold off;

            pause(2);
            exportgraphics(gca, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);
            
    end
    
end

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ')
disp('***********');

disp(' ');


%% Clean

clearvars -except pressureData
disp(' ');