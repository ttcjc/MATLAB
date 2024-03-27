run preamble;

load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_120s_15Hz_02/T0067_T120000_F15');
% load('/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/Far_Field_Soiling_07_22/SB_1.0L_600s_03Hz_01/T0333_T600000_F3.mat');

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;

xLimsPlot = [0.3; 4.6257662];
yLimsPlot = [-0.25; 0.25];
zLimsPlot = [0; 0.4];


%%

yLims = [-0.2; 0.1];
zLims = [min(mapData.positionGrid(:,3) / normLength); zDims(2)];

y = yLims(1):(diff(yLims) / 10):yLims(2);
z = zLims(1):(diff(zLims) / 10):zLims(2);

[y, z] = ndgrid(y, z);

samplePoints = [y(:), z(:)];

index = [12; 17; 22; 61; 116];

samplePoints = samplePoints(index,:);

index = dsearchn((mapData.positionGrid(:,[2,3]) / normLength), samplePoints);


%%

nTimes = height(mapData.time);

%%%

expandingMean = cell(nTimes, 1); expandingMean(:) = {zeros([height(index), 1])};

for i = 1:nTimes
    
    if i == 1
        expandingMean{i} = full(mapData.density.inst{i}(index));
    else
        expandingMean{i} = expandingMean{i-1} + full(mapData.density.inst{i}(index));
    end
    
end

for i = 1:nTimes
    expandingMean{i} = expandingMean{i} / i;
end

%%%

% error = cell(nTimes, 1); error(:) = {zeros([1, height(index)])};
% 
% for i = 1:nTimes
%     error{i} = expandingMean{i}';
% end
% 
% error = cell2mat(error);

%%%

error = cell((nTimes - 1), 1); error(:) = {zeros([1, height(index)])};

for i = 2:nTimes
    error{i - 1} = abs(expandingMean{i} - expandingMean{i - 1})';
end

error = cell2mat(error);


%%

% Initialise Figure
fig = fig + 1;
figName = 'Convergence_of_Time_Average';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', 'yScale', 'log');

% Plot Mass Flux
plot((movmean(error(:,1), 8) / mean(mapData.density.mean(mapData.density.mean > 0))), ...
     'color', graphColours(1), 'lineWidth', 2);
plot((movmean(error(:,2), 8) / mean(mapData.density.mean(mapData.density.mean > 0))), ...
     'color', graphColours(2), 'lineWidth', 2);
plot((movmean(error(:,3), 8) / mean(mapData.density.mean(mapData.density.mean > 0))), ...
     'color', graphColours(3), 'lineWidth', 2);
plot((movmean(error(:,4), 8) / mean(mapData.density.mean(mapData.density.mean > 0))), ...
     'color', graphColours(4), 'lineWidth', 2);
plot((movmean(error(:,5), 8) / mean(mapData.density.mean(mapData.density.mean > 0))), ...
     'color', graphColours(5), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 1800]);
ylim([1e-5; 1e0]);
tickData = (360:360:1440);
xticks(tickData);
tickData = [1e-4; 1e-3; 1e-2; 1e-1];
yticks(tickData);
xlabel({'{Samples}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\epsilon$}'}, 'interpreter', 'latex');
legend({'Point A', ...
        'Point B', ...
        'Point C', ...
        'Point D', ...
        'Point E'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');