clA;

content = importdata('~/Data/Particle Size Data/Bouchet Data/80 kph.csv');


%% 

[diameterRaw, index] = unique(content(:,1)); diameterRaw = diameterRaw * 1e3;
populationRaw = content(index,2);

%%

minD = 20; % um
maxD = 400; % um

interpSize = 1; % um
diameterInterp = ((minD - (interpSize / 2)):(interpSize / 2):(maxD + (interpSize / 2)))';
populationInterp = interp1(diameterRaw, populationRaw, diameterInterp, 'linear', 'extrap');


%%

interpolant = griddedInterpolant(diameterInterp, populationInterp);

binSize = 2; % um

diameterBinned = (minD:binSize:maxD)';
populationBinned = zeros([height(diameterBinned),1]);

for i = 1:height(diameterBinned)
    a = interpolant(diameterBinned(i) - (binSize / 2));
    b = interpolant(diameterBinned(i) - (binSize / 2));

    populationBinned(i) = ((a + b) / 2) * binSize;
end

populationBinned = (populationBinned / sum(populationBinned)) * 100;


%%

% Initialise Figure
set(figure(1), 'name', 'Bouchet Distribution', 'color', [1, 1, 1], ...
                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top');
hold on;

% Plot Data
plot(diameterInterp, populationInterp, 'color', 'r', 'lineStyle', '-', 'lineWidth', 2);
scatter(diameterRaw, populationRaw, 50, 'r');
scatter(diameterBinned, populationBinned, 50, 'b');
hold off;
xlim([0, 420]);
ylim([0, 10]);


%%

fileID = fopen('~/MATLAB/Output/Files/Bouchet_Distribution', 'w');
formatSpec = '                        (%e    %e)\n';

for i = 1:height(diameterBinned)
    fprintf(fileID, formatSpec, (diameterBinned(i) * 1e-6), populationBinned(i));
end

fclose(fileID);