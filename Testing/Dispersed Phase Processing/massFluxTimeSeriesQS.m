run preamble;

SB = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawQS_SB.mat', 'LagData').LagData;
ST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawQS_ST.mat', 'LagData').LagData;
RSST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawQS_RSST.mat', 'LagData').LagData;

planes = fieldnames(SB);

plotPlane = 4;

dT = 2e-5;


%% Calculate Instantaneous Mass Flux

for i = 1:height(planes)
    nTimes = height(SB.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = SB.(planes{i}).nParticle;
    d = SB.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    SB.(planes{i}).massFlow = massFlow; clear massFlow;
    
    %%%
    
    nTimes = height(ST.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = ST.(planes{i}).nParticle;
    d = ST.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    ST.(planes{i}).massFlow = massFlow; clear massFlow;
    
    %%%
    
    nTimes = height(RSST.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = RSST.(planes{i}).nParticle;
    d = RSST.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    RSST.(planes{i}).massFlow = massFlow; clear massFlow;
end
clear i;


%% Plot Instantaneous Data

% Initialise Figure
fig = fig + 1;
figName = ['Instantaneous_Mass_Flux_', num2str(plotPlane)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Instantaneous Mass Flux
plot(SB.(planes{plotPlane}).time, (SB.(planes{plotPlane}).massFlow / dT), ...
     'color', graphColours(1), 'lineWidth', 2);

plot(SB.(planes{plotPlane}).time, movmean((SB.(planes{plotPlane}).massFlow / dT), 512), ...
     'color', graphColours(2), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1; 4]);
ylim([0; 3.5e-3]);
tickData = (1.6:0.6:3.4);
xticks(tickData);
tickData = (0.7e-3:0.7e-3:2.8e-3);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.1f');
currentFig = gca; currentFig.YAxis.Exponent = -3; clear currentAxis;
xlabel({'{Time $(s)$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\dot{m}_{_{p}}$ $(kg {\cdot} s^{-1})$}'}, 'interpreter', 'latex');
legend({'Instantaneous', ...
        'Moving Mean'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Plot Cumulative Data

% Initialise Figure
fig = fig + 1;
figName = 'Cumulative_Mass_Flux';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Cumulative Mass Flux
plot(SB.(planes{plotPlane}).time, cumsum(SB.(planes{plotPlane}).massFlow), ...
     'color', graphColours(1), 'lineWidth', 2);
plot(ST.(planes{plotPlane}).time, cumsum(ST.(planes{plotPlane}).massFlow), ...
     'color', graphColours(2), 'lineWidth', 2);
plot(RSST.(planes{plotPlane}).time, cumsum(RSST.(planes{plotPlane}).massFlow), ...
     'color', graphColours(3), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1; 4]);
ylim([0; 3.5e-3]);
tickData = (1.6:0.6:3.4);
xticks(tickData);
tickData = (0.7e-3:0.7e-3:2.8e-3);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.1f');
currentFig = gca; currentFig.YAxis.Exponent = -3; clear currentAxis;
xlabel({'{Time $(s)$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Cumulative Mass Transfer $(kg)$}'}, 'interpreter', 'latex');
legend({'\textit{Config A}', ...
        '\textit{Config B}', ...
        '\textit{Config C}'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Plot Total Data

totalMassSB = zeros([height(planes),1], 'single');
totalMassST = totalMassSB;
totalMassRSST = totalMassSB;

for i = 1:height(planes)
    
    if i == 1
        totalMassSB(i) = NaN;
        totalMassST(i) = NaN;
        totalMassRSST(i) = NaN;
    else
        totalMassSB(i) = max(cumsum(SB.(planes{i}).massFlow));
        totalMassST(i) = max(cumsum(ST.(planes{i}).massFlow));
        totalMassRSST(i) = max(cumsum(RSST.(planes{i}).massFlow));
    end
    
end

% Initialise Figure
fig = fig + 1;
figName = 'Total_Mass_Flux';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');


% Plot Mass
plot(totalMassSB,   'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2, ...
                    'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(1));
plot(totalMassST,   'color', graphColours(2), 'lineStyle', '-', 'lineWidth', 2, ...
                    'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(2));
plot(totalMassRSST, 'color', graphColours(3), 'lineStyle', '-', 'lineWidth', 2, ...
                    'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(3));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1.5; 4.5]);
ylim([3e-3; 4e-3]);
tickData = [2; 3; 4];
xticks(tickData);
tickData = (3.2e-3:0.2e-3:3.8e-3);
yticks(tickData);
currentAxisX = get(gca, 'XAxis'); currentAxisX.TickLabelInterpreter = 'latex'; clear currentAxisX;
set(gca, 'XTickLabel', {'$1.0\,\ell$', '$1.5\,\ell$', '$2.0\,\ell$'});
ytickformat('%.1f');
currentFig = gca; currentFig.YAxis.Exponent = -3; clear currentAxis;
xlabel({'{Measurement Plane}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Total Mass Flow ($kg$)}'}, 'interpreter', 'latex');
legend({'\textit{Config A}', ...
        '\textit{Config B}', ...
        '\textit{Config C}'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
