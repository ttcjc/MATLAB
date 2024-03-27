run preamble;

% SB = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawSB.mat', 'LagData').LagData;
% ST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawST.mat', 'LagData').LagData;
% RSST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawRSST.mat', 'LagData').LagData;
% 
% planes = fieldnames(SB);


uncoupled = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_uncoupled.mat', 'LagData').LagData;
coupled = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_coupled.mat', 'LagData').LagData;
halfTread = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_halfTread.mat', 'LagData').LagData;
twentyDeg = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_20deg.mat', 'LagData').LagData;

uncoupled = orderfields(uncoupled, [4,1,2,3]);
coupled = orderfields(coupled, [4,1,2,3]);
halfTread = orderfields(halfTread, [4,1,2,3]);
twentyDeg = orderfields(twentyDeg, [4,1,2,3]);

planes = fieldnames(uncoupled);

plotPlane = 1;


%% Calculate Instantaneous Mass Flux

for i = 1:height(planes)
    nTimes = height(uncoupled.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = uncoupled.(planes{i}).nParticle;
    d = uncoupled.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    uncoupled.(planes{i}).massFlow = massFlow; clear massFlow;
    
    %%%
    
    nTimes = height(coupled.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = coupled.(planes{i}).nParticle;
    d = coupled.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    coupled.(planes{i}).massFlow = massFlow; clear massFlow;
    
    %%%
    
    nTimes = height(halfTread.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = halfTread.(planes{i}).nParticle;
    d = halfTread.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    halfTread.(planes{i}).massFlow = massFlow; clear massFlow;
    
    %%%
    
    nTimes = height(twentyDeg.(planes{i}).time);
    
    massFlow = zeros([nTimes, 1], 'single');

    nParticle = twentyDeg.(planes{i}).nParticle;
    d = twentyDeg.(planes{i}).d;
    for j = 1:nTimes
        massFlow(j) = sum(nParticle{j} .* (1000 * ((1 / 12) * tau * d{j}.^3)));
    end
    clear j nParticle d;

    twentyDeg.(planes{i}).massFlow = massFlow; clear massFlow;
end
clear i;


%% Plot Instantaneous Data

% % Initialise Figure
% fig = fig + 1;
% figName = 'Instantaneous_Mass_Flux';
% set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%              'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
% pause(0.5);
% hold on;
% set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%          'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
% 
% % Plot Instantaneous Mass Flux
% plot(SB.(planes{end}).time, SB.(planes{end}).massFlow, ...
%      'color', graphColours(1), 'lineWidth', 2);
% plot(ST.(planes{end}).time, ST.(planes{end}).massFlow, ...
%      'color', graphColours(2), 'lineWidth', 2);
% plot(RSST.(planes{end}).time, RSST.(planes{end}).massFlow, ...
%      'color', graphColours(3), 'lineWidth', 2);
% 
% % Format Figure
% title('{-----}', 'interpreter', 'latex');
% subtitle('{ }');
% axis on;
% box on;
% grid off;
% xlim([1; 4]);
% ylim([0; 70e-9]);
% tickData = (1.6:0.6:3.4);
% xticks(tickData);
% tickData = (14e-9:14e-9:56e-9);
% yticks(tickData);
% % xtickformat('%.1f');
% % ytickformat('%.1f');
% xlabel({'{Time ($s$)}'; '{-----}'}, 'interpreter', 'latex');
% ylabel({'{-----}'; '{Mass ($kg$)}'}, 'interpreter', 'latex');
% legend({'\textit{Config A}', ...
%         '\textit{Config B}', ...
%         '\textit{Config C}'}, ...
%        'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
%        'fontSize', 18, 'box', 'off');
% tightInset = get(gca, 'TightInset');
% set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                            (tightInset(2) + 0.00625), ...
%                            (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                            (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
% pause(0.5);
% hold off;
% 
% % Save Figure
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');

% Initialise Figure
fig = fig + 1;
figName = ['Instantaneous_Mass_Flux_', num2str(plotPlane)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% % Plot Instantaneous Mass Flux
% plot(uncoupled.(planes{plotPlane}).time, uncoupled.(planes{plotPlane}).massFlow, ...
%      'color', graphColours(4), 'lineWidth', 2);
% plot(coupled.(planes{plotPlane}).time, coupled.(planes{plotPlane}).massFlow, ...
%      'color', graphColours(5), 'lineWidth', 2);
% plot(halfTread.(planes{plotPlane}).time, halfTread.(planes{plotPlane}).massFlow, ...
%      'color', graphColours(6), 'lineWidth', 2);
% plot(twentyDeg.(planes{plotPlane}).time, twentyDeg.(planes{plotPlane}).massFlow, ...
%      'color', graphColours(7), 'lineWidth', 2);

% Plot Moving Mean
plot(uncoupled.(planes{plotPlane}).time, movmean(uncoupled.(planes{plotPlane}).massFlow, 1048), ...
     'color', graphColours(4), 'lineWidth', 2);
plot(coupled.(planes{plotPlane}).time, movmean(coupled.(planes{plotPlane}).massFlow, 1048), ...
     'color', graphColours(5), 'lineWidth', 2);
plot(halfTread.(planes{plotPlane}).time, movmean(halfTread.(planes{plotPlane}).massFlow, 1048), ...
     'color', graphColours(6), 'lineWidth', 2);
plot(twentyDeg.(planes{plotPlane}).time, movmean(twentyDeg.(planes{plotPlane}).massFlow, 1048), ...
     'color', graphColours(7), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([10; 32]);
ylim([0; 3e-3]);
tickData = (14.4:4.4:27.6);
xticks(tickData);
tickData = (0.6e-3:0.6e-3:2.4e-3);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.1f');
xlabel({'{Time ($s$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Mass ($kg$)}'}, 'interpreter', 'latex');
legend({'Uncoupled', ...
        'Coupled', ...
        'Reduced Mass', ...
        'Increased Angle'}, ...
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
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Plot Cumulative Data

% % Initialise Figure
% fig = fig + 1;
% figName = 'Cumulative_Mass_Flux';
% set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%              'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
% pause(0.5);
% hold on;
% set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%          'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
% 
% % Plot Cumulative Mass Flux
% plot(SB.(planes{end}).time, cumsum(SB.(planes{end}).massFlow), ...
%      'color', graphColours(1), 'lineWidth', 2);
% plot(ST.(planes{end}).time, cumsum(ST.(planes{end}).massFlow), ...
%      'color', graphColours(2), 'lineWidth', 2);
% plot(RSST.(planes{end}).time, cumsum(RSST.(planes{end}).massFlow), ...
%      'color', graphColours(3), 'lineWidth', 2);
% 
% % Format Figure
% title('{-----}', 'interpreter', 'latex');
% subtitle('{ }');
% axis on;
% box on;
% grid off;
% xlim([1; 4]);
% ylim([0; 3.5e-3]);
% tickData = (1.6:0.6:3.4);
% xticks(tickData);
% tickData = (0.7e-3:0.7e-3:2.8e-3);
% yticks(tickData);
% xtickformat('%.1f');
% ytickformat('%.1f');
% xlabel({'{Time ($s$)}'; '{-----}'}, 'interpreter', 'latex');
% ylabel({'{-----}'; '{Mass ($kg$)}'}, 'interpreter', 'latex');
% legend({'\textit{Config A}', ...
%         '\textit{Config B}', ...
%         '\textit{Config C}'}, ...
%        'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
%        'fontSize', 18, 'box', 'off');
% tightInset = get(gca, 'TightInset');
% set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                            (tightInset(2) + 0.00625), ...
%                            (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                            (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
% pause(0.5);
% hold off;
% 
% % Save Figure
% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


% Initialise Figure
fig = fig + 1;
figName = ['Cumulative_Mass_Flux_', num2str(plotPlane)];
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Cumulative Mass Flux
plot(uncoupled.(planes{plotPlane}).time, cumsum(uncoupled.(planes{plotPlane}).massFlow), ...
     'color', graphColours(4), 'lineWidth', 2);
plot(coupled.(planes{plotPlane}).time, cumsum(coupled.(planes{plotPlane}).massFlow), ...
     'color', graphColours(5), 'lineWidth', 2);
plot(halfTread.(planes{plotPlane}).time, cumsum(halfTread.(planes{plotPlane}).massFlow), ...
     'color', graphColours(6), 'lineWidth', 2);
plot(twentyDeg.(planes{plotPlane}).time, cumsum(twentyDeg.(planes{plotPlane}).massFlow), ...
     'color', graphColours(7), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([10; 32]);
ylim([0; 250]);
tickData = (14.4:4.4:27.6);
xticks(tickData);
tickData = (50:50:200);
yticks(tickData);
xtickformat('%.1f');
ytickformat('%.0f');
xlabel({'{Time ($s$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Mass ($kg$)}'}, 'interpreter', 'latex');
legend({'Uncoupled', ...
        'Coupled', ...
        'Reduced Mass', ...
        'Increased Angle'}, ...
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
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%%

% startIndex = find(uncoupled.(planes{plotPlane}).time == single(10.000125));
% 
% nTimes = height(uncoupled.(planes{plotPlane}).time(startIndex:end));
% 
% massMean_Uncoupled = mean(uncoupled.(planes{plotPlane}).massFlow(startIndex:end));
% 
% massPrime_Uncoupled = uncoupled.(planes{plotPlane}).massFlow(startIndex:end) - massMean_Uncoupled;
% 
% massPrimeNorm_Uncoupled = massPrime_Uncoupled ./ massMean_Uncoupled;
% 
% massFlowMean_Uncoupled = massMean_Uncoupled / 1.25e-4
% 
% massRMS_Uncoupled = sqrt((1 / nTimes) * sum(massPrimeNorm_Uncoupled.^2))
% 
% massDeviation_Uncoupled = prctile(abs(massPrimeNorm_Uncoupled), 99)
% 
% 
% %%%
% 
% startIndex = find(coupled.(planes{plotPlane}).time == single(10.000125));
% 
% nTimes = height(coupled.(planes{plotPlane}).time(startIndex:end));
% 
% massMean_Coupled = mean(coupled.(planes{plotPlane}).massFlow(startIndex:end));
% 
% massPrime_Coupled = coupled.(planes{plotPlane}).massFlow(startIndex:end) - massMean_Coupled;
% 
% massPrimeNorm_Coupled = massPrime_Coupled ./ massMean_Coupled;
% 
% massFlowMean_Coupled = massMean_Coupled / 1.25e-4
% 
% massRMS_Coupled = sqrt((1 / nTimes) * sum(massPrimeNorm_Coupled.^2))
% 
% massDeviation_Coupled = prctile(abs(massPrimeNorm_Coupled), 99)


%%

totalMassUncoupled = zeros([height(planes),1], 'single');
totalMassHalfTread = totalMassUncoupled;
totalMassCoupled = totalMassUncoupled;
totalMass20deg = totalMassUncoupled;

for i = 1:height(planes)
    totalMassUncoupled(i) = max(cumsum(uncoupled.(planes{i}).massFlow));
    totalMassHalfTread(i) = max(cumsum(halfTread.(planes{i}).massFlow));
    totalMassCoupled(i) = max(cumsum(coupled.(planes{i}).massFlow));
    totalMass20deg(i) = max(cumsum(twentyDeg.(planes{i}).massFlow));
end

% Initialise Figure
fig = fig + 1;
figName = 'Total_Mass_Through_Planes';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');


% Plot Mass
plot(totalMassUncoupled, 'color', graphColours(4), 'lineStyle', '-', 'lineWidth', 2, ...
                         'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(4));
plot(totalMassCoupled,   'color', graphColours(5), 'lineStyle', '-', 'lineWidth', 2, ...
                         'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(5));
plot(totalMassHalfTread, 'color', graphColours(6), 'lineStyle', '-', 'lineWidth', 2, ...
                         'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(6));
plot(totalMass20deg,     'color', graphColours(7), 'lineStyle', '-', 'lineWidth', 2, ...
                         'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(7));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 5]);
ylim([0; 250]);
tickData = (1:1:4);
xticks(tickData);
set(gca, 'XTickLabel', {'$1\,\ell$', '$2\,\ell$', '$3\,\ell$', '$4\,\ell$'}, 'TickLabelInterpreter', 'latex');
tickData = (50:50:200);
yticks(tickData);
xlabel({'{Measurement Plane}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Total Mass ($kg$)}'}, 'interpreter', 'latex');
legend({'Uncoupled', ...
        'Coupled', ...
        'Reduced Mass', ...
        'Increased Angle'}, ...
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