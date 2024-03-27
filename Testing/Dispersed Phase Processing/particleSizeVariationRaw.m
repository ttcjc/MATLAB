run preamble;

SB = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawSB.mat', 'LagData').LagData;
ST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawST.mat', 'LagData').LagData;
RSST = load('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawRSST.mat', 'LagData').LagData;

planes = fieldnames(SB);


%%

dDist = single((1e-6:1e-6:147e-6)');

sizeDist.injector = [
                    1.000000e-06, 1.499701e-02
                    2.000000e-06, 8.915905e-02
                    3.000000e-06, 1.265001e-01
                    4.000000e-06, 1.510801e-01
                    5.000000e-06, 1.895401e-01
                    6.000000e-06, 2.475902e-01
                    7.000000e-06, 3.372802e-01
                    8.000000e-06, 4.288903e-01
                    9.000000e-06, 5.191703e-01
                    1.000000e-05, 6.351104e-01
                    1.100000e-05, 7.801405e-01
                    1.200000e-05, 9.656206e-01
                    1.300000e-05, 1.110601e+00
                    1.400000e-05, 1.340601e+00
                    1.500000e-05, 1.644301e+00
                    1.600000e-05, 1.934001e+00
                    1.700000e-05, 2.221901e+00
                    1.800000e-05, 2.465402e+00
                    1.900000e-05, 2.700402e+00
                    2.000000e-05, 2.956402e+00
                    2.100000e-05, 3.169502e+00
                    2.200000e-05, 3.328702e+00
                    2.300000e-05, 3.473002e+00
                    2.400000e-05, 3.557502e+00
                    2.500000e-05, 3.642902e+00
                    2.600000e-05, 3.611302e+00
                    2.700000e-05, 3.665702e+00
                    2.800000e-05, 3.586202e+00
                    2.900000e-05, 3.538302e+00
                    3.000000e-05, 3.468802e+00
                    3.100000e-05, 3.308302e+00
                    3.200000e-05, 3.150302e+00
                    3.300000e-05, 3.061802e+00
                    3.400000e-05, 2.829402e+00
                    3.500000e-05, 2.640902e+00
                    3.600000e-05, 2.542902e+00
                    3.700000e-05, 2.353501e+00
                    3.800000e-05, 2.221601e+00
                    3.900000e-05, 1.983301e+00
                    4.000000e-05, 1.756801e+00
                    4.100000e-05, 1.694401e+00
                    4.200000e-05, 1.565901e+00
                    4.300000e-05, 1.369001e+00
                    4.400000e-05, 1.291501e+00
                    4.500000e-05, 1.189001e+00
                    4.600000e-05, 1.063601e+00
                    4.700000e-05, 9.799506e-01
                    4.800000e-05, 8.573805e-01
                    4.900000e-05, 7.601505e-01
                    5.000000e-05, 7.387305e-01
                    5.100000e-05, 6.354004e-01
                    5.200000e-05, 5.984604e-01
                    5.300000e-05, 5.149103e-01
                    5.400000e-05, 4.849603e-01
                    5.500000e-05, 4.227703e-01
                    5.600000e-05, 3.829702e-01
                    5.700000e-05, 3.633602e-01
                    5.800000e-05, 2.899102e-01
                    5.900000e-05, 2.793302e-01
                    6.000000e-05, 2.576302e-01
                    6.100000e-05, 2.447202e-01
                    6.200000e-05, 2.166601e-01
                    6.300000e-05, 1.917901e-01
                    6.400000e-05, 1.537601e-01
                    6.500000e-05, 1.587401e-01
                    6.600000e-05, 1.139901e-01
                    6.700000e-05, 1.230401e-01
                    6.800000e-05, 1.152001e-01
                    6.900000e-05, 9.778006e-02
                    7.000000e-05, 8.969306e-02
                    7.100000e-05, 7.103304e-02
                    7.200000e-05, 7.348205e-02
                    7.300000e-05, 6.799404e-02
                    7.400000e-05, 6.054504e-02
                    7.500000e-05, 7.010404e-02
                    7.600000e-05, 5.180003e-02
                    7.700000e-05, 4.161203e-02
                    7.800000e-05, 3.733802e-02
                    7.900000e-05, 5.876904e-02
                    8.000000e-05, 3.205302e-02
                    8.100000e-05, 2.832802e-02
                    8.200000e-05, 3.885702e-02
                    8.300000e-05, 2.323901e-02
                    8.400000e-05, 2.494202e-02
                    8.500000e-05, 3.339602e-02
                    8.600000e-05, 3.349302e-02
                    8.700000e-05, 1.232601e-02
                    8.800000e-05, 1.750001e-02
                    8.900000e-05, 1.255001e-02
                    9.000000e-05, 1.230101e-02
                    9.100000e-05, 1.479801e-02
                    9.200000e-05, 1.157701e-02
                    9.300000e-05, 1.134501e-02
                    9.400000e-05, 6.837504e-03
                    9.500000e-05, 1.011101e-02
                    9.600000e-05, 1.012001e-02
                    9.700000e-05, 9.452306e-03
                    9.800000e-05, 1.119401e-02
                    9.900000e-05, 1.216401e-02
                    1.000000e-04, 9.294306e-03
                    1.010000e-04, 7.158404e-03
                    1.020000e-04, 3.451902e-03
                    1.030000e-04, 7.978905e-03
                    1.040000e-04, 1.359201e-02
                    1.050000e-04, 1.047601e-02
                    1.060000e-04, 7.137904e-03
                    1.070000e-04, 2.752102e-03
                    1.080000e-04, 1.731201e-03
                    1.090000e-04, 1.796701e-03
                    1.100000e-04, 4.850803e-03
                    1.120000e-04, 1.513101e-03
                    1.130000e-04, 5.535903e-03
                    1.140000e-04, 2.893402e-03
                    1.170000e-04, 2.563102e-03
                    1.180000e-04, 4.448903e-03
                    1.190000e-04, 3.210302e-03
                    1.210000e-04, 1.077601e-03
                    1.250000e-04, 1.326001e-03
                    1.270000e-04, 1.357701e-03
                    1.290000e-04, 4.128703e-03
                    1.320000e-04, 1.686201e-03
                    1.330000e-04, 2.157601e-03
                    1.340000e-04, 2.217001e-04
                    1.350000e-04, 1.941201e-04
                    1.360000e-04, 2.426901e-04
                    1.370000e-04, 2.162801e-04
                    1.390000e-04, 4.796303e-04
                    1.400000e-04, 6.741304e-04
                    1.410000e-04, 1.820201e-04
                    1.420000e-04, 8.302505e-04
                    1.430000e-04, 1.491001e-04
                    1.440000e-04, 7.650305e-04
                    1.460000e-04, 1.694801e-04
                    1.470000e-04, 3.754902e-04
                   ];


%% SB

nTimes = height(SB.(planes{end}).time);

% Calculate Instantaneous Distributions
sizeDist.SB.(planes{end}).inst = cell(nTimes, 1); sizeDist.SB.(planes{end}).inst(:) = {[dDist, zeros([height(dDist), 1], 'single')]};

for j = 1:nTimes
    massParcel = SB.(planes{end}).nParticle{j} .* ...
                 (1000 * ((1 / 12) * tau * SB.(planes{end}).d{j}.^3));

    dBinned = round(SB.(planes{end}).d{j}, 6);

    for k = 1:height(dDist)
        index = find(dBinned == dDist(k));

        if ~isempty(index)
            sizeDist.SB.(planes{end}).inst{j}(k,2) = sum(massParcel(index)) / sum(massParcel);
        end

    end
    clear k;

end
clear j;

% Calculate Time-Averaged Distributions
sizeDist.SB.(planes{end}).mean = [dDist, zeros([height(dDist), 1], 'single')];

for j = 1:nTimes
    sizeDist.SB.(planes{end}).mean(:,2) = sizeDist.SB.(planes{end}).mean(:,2) + sizeDist.SB.(planes{end}).inst{j}(:,2);
end
clear j;

sizeDist.SB.(planes{end}).mean((sizeDist.SB.(planes{end}).mean(:,2) == 0), 2) = NaN;

sizeDist.SB.(planes{end}).mean(:,2) = sizeDist.SB.(planes{end}).mean(:,2) / nTimes;


%% ST

nTimes = height(ST.(planes{end}).time);

% Calculate Instantaneous Distributions
sizeDist.ST.(planes{end}).inst = cell(nTimes, 1); sizeDist.ST.(planes{end}).inst(:) = {[dDist, zeros([height(dDist), 1], 'single')]};

for j = 1:nTimes
    massParcel = ST.(planes{end}).nParticle{j} .* ...
                 (1000 * ((1 / 12) * tau * ST.(planes{end}).d{j}.^3));

    dBinned = round(ST.(planes{end}).d{j}, 6);

    for k = 1:height(dDist)
        index = find(dBinned == dDist(k));

        if ~isempty(index)
            sizeDist.ST.(planes{end}).inst{j}(k,2) = sum(massParcel(index)) / sum(massParcel);
        end

    end
    clear k;

end
clear j;

% Calculate Time-Averaged Distributions
sizeDist.ST.(planes{end}).mean = [dDist, zeros([height(dDist), 1], 'single')];

for j = 1:nTimes
    sizeDist.ST.(planes{end}).mean(:,2) = sizeDist.ST.(planes{end}).mean(:,2) + sizeDist.ST.(planes{end}).inst{j}(:,2);
end
clear j;

sizeDist.ST.(planes{end}).mean((sizeDist.ST.(planes{end}).mean(:,2) == 0), 2) = NaN;

sizeDist.ST.(planes{end}).mean(:,2) = sizeDist.ST.(planes{end}).mean(:,2) / nTimes;


%% RSST

nTimes = height(RSST.(planes{end}).time);

% Calculate Instantaneous Distributions
sizeDist.RSST.(planes{end}).inst = cell(nTimes, 1); sizeDist.RSST.(planes{end}).inst(:) = {[dDist, zeros([height(dDist), 1], 'single')]};

for j = 1:nTimes
    massParcel = RSST.(planes{end}).nParticle{j} .* ...
                 (1000 * ((1 / 12) * tau * RSST.(planes{end}).d{j}.^3));

    dBinned = round(RSST.(planes{end}).d{j}, 6);

    for k = 1:height(dDist)
        index = find(dBinned == dDist(k));

        if ~isempty(index)
            sizeDist.RSST.(planes{end}).inst{j}(k,2) = sum(massParcel(index)) / sum(massParcel);
        end

    end
    clear k;

end
clear j;

% Calculate Time-Averaged Distributions
sizeDist.RSST.(planes{end}).mean = [dDist, zeros([height(dDist), 1], 'single')];

for j = 1:nTimes
    sizeDist.RSST.(planes{end}).mean(:,2) = sizeDist.RSST.(planes{end}).mean(:,2) + sizeDist.RSST.(planes{end}).inst{j}(:,2);
end
clear j;

sizeDist.RSST.(planes{end}).mean((sizeDist.RSST.(planes{end}).mean(:,2) == 0), 2) = NaN;

sizeDist.RSST.(planes{end}).mean(:,2) = sizeDist.RSST.(planes{end}).mean(:,2) / nTimes;


%% Plot Data

% Initialise Figure
fig = fig + 1;
figName = 'Particle_Size_Variation';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', 'yScale', 'log');

% Plot Particle Size Distributions
plot((sizeDist.injector(:,1) * 1e6), sizeDist.injector(:,2), ...
     'color', graphColours(4), 'lineWidth', 2);      
plot((sizeDist.SB.(planes{end}).mean(:,1) * 1e6), (sizeDist.SB.(planes{end}).mean(:,2) * 100), ...
     'color', graphColours(1), 'lineWidth', 2);
plot((sizeDist.ST.(planes{end}).mean(:,1) * 1e6), (sizeDist.ST.(planes{end}).mean(:,2) * 100), ...
     'color', graphColours(2), 'lineWidth', 2);
plot((sizeDist.RSST.(planes{end}).mean(:,1) * 1e6), (sizeDist.RSST.(planes{end}).mean(:,2) * 100), ...
     'color', graphColours(3), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 150]);
ylim([1e-4; 1e1]);
tickData = (30:30:120);
xticks(tickData);
tickData = [1e-3; 1e-2; 1e-1; 1e0];
yticks(tickData);
xtickformat('%.0f');
xlabel({'{$D_{_{p}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Mass-Weighted Population (\%)}'}, 'interpreter', 'latex');
% legend({'Injector', ...
legend({'Injector', ...
        '\textit{Config A}', ...
        '\textit{Config B}', ...
        '\textit{Config C}'}, ...
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