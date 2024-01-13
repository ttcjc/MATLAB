%% Lagrangian Mass Flux Time Series
% ----
% Lorem ipsum


%% Preamble

run preamble;

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

normMass = true; % Normalise Mass in Plots
    normValue = 0.0034302; % Windsor_SB_wW_Upstream_SC 1.0L

disp('================================');
disp('Lagrangian Mass Flux Time Series');
disp('================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Case

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData(geoLoc);

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

maxNumCompThreads(nProc);

[dataID, LagProps, ~, LagData, ...
 ~, sampleInt, dataFormat] = initialiseLagData(saveLoc, caseFolder, campaignID, ...
                                               caseID, cloudName, false, true, ...
                                               false, timeDirs, deltaT, ...
                                               timePrecision, maxNumCompThreads);

nTimes = height(LagData.X_N0_31975.time);

planes = fieldnames(LagData);


%% Calculate Planar Mass Flux

for i = 2:height(planes)
    massFlux.(planes{i}).inst = zeros([nTimes, 1], 'single');
    
    for j = 1:nTimes
        massFlux.(planes{i}).inst(j) = sum(LagData.(planes{i}).nParticle{j} .* ...
                                           (1000 .* ((1 / 12) * tau * LagData.(planes{i}).d{j}.^3)));
    end
    clear j;
    
end
clear i;


%% Plot Data

% Initialise Figure
fig = fig + 1;
figName = ['Cumulative Mass Flux_', caseID];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mass Flux
for i = 2:height(planes)
    plot(LagData.(planes{2}).time, (cumsum(massFlux.(planes{i}).inst) / normValue), ...
         'color', graphColours(i - 1), 'lineWidth', 2);   
end
clear i;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1; 4]);
ylim([0; 1.05]);
tickData = (1.6:0.6:3.4);
xticks(tickData);
tickData = (0.21:0.21:0.84);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.2f');
xlabel({'{Time ($s$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$m_{_{n}}$}'}, 'interpreter', 'latex');
legend({'1.0\,$\ell$ Measurement Plane', '1.5\,$\ell$ Measurement Plane', '2.0\,$\ell$ Measurement Plane'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off')
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');