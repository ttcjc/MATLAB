run preamble;

load('~/MATLAB/Testing/coeffData_QS_SB.mat');

startTime = 1;
startSample = find(coeffData.time == startTime);

coeffData.time = coeffData.time(startSample:end);


%%

% Perform FFT
coeffData.Cd = coeffData.Cd(startSample:end);

n = pow2(nextpow2(height(coeffData.time)));
Fs = round(1 / (coeffData.time(2) - coeffData.time(1)));

x = zeros(n,1);
x(1:height(coeffData.Cd)) = coeffData.Cd - mean(coeffData.Cd);

[PSD, freq] = pwelch(x, (n / 8), [], n, Fs);
Sr = (freq * 1.044) / 40;

% Figure Setup
fig = fig + 1;
figName = 'Drag_Coefficient_Temporal_Analysis';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', ...
         'xScale', 'log', 'yScale', 'log');

% Plot FFT
plot(Sr(1:floor(n / 2)), rescale(PSD(1:floor(n / 2))), ...
     'color', graphColours(1), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1e-2; 1e3]);
ylim([1e-13; 1e2]);
tickData = [1e-1; 1e0; 1e1; 1e2];
xticks(tickData);
tickData = [1e-10; 1e-7; 1e-4; 1e-1];
yticks(tickData);
% xtickformat('%.1f');
% ytickformat('%.1f');
xlabel({'{$St_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{PSD}'}, 'interpreter', 'latex');
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

% Perform FFT
coeffData.Cl = coeffData.Cl(startSample:end);

n = pow2(nextpow2(height(coeffData.time)));
Fs = round(1 / (coeffData.time(2) - coeffData.time(1)));

x = zeros(n,1);
x(1:height(coeffData.Cl)) = coeffData.Cl - mean(coeffData.Cl);

[PSD, freq] = pwelch(x, (n / 8), [], n, Fs);
Sr = (freq * 1.044) / 40;

% Figure Setup
fig = fig + 1;
figName = 'Lift_Coefficient_Temporal_Analysis';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', ...
         'xScale', 'log', 'yScale', 'log');

% Plot FFT
plot(Sr(1:floor(n / 2)), rescale(PSD(1:floor(n / 2))), ...
     'color', graphColours(1), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1e-2; 1e3]);
ylim([1e-13; 1e2]);
tickData = [1e-1; 1e0; 1e1; 1e2];
xticks(tickData);
tickData = [1e-10; 1e-7; 1e-4; 1e-1];
yticks(tickData);
% xtickformat('%.1f');
% ytickformat('%.1f');
xlabel({'{$St_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{PSD}'}, 'interpreter', 'latex');
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

% Perform FFT
coeffData.Cs = coeffData.Cs(startSample:end);

n = pow2(nextpow2(height(coeffData.time)));
Fs = round(1 / (coeffData.time(2) - coeffData.time(1)));

x = zeros(n,1);
x(1:height(coeffData.Cs)) = coeffData.Cs - mean(coeffData.Cs);

[PSD, freq] = pwelch(x, (n / 8), [], n, Fs);
Sr = (freq * 1.044) / 40;

% Figure Setup
fig = fig + 1;
figName = 'Side_Coefficient_Temporal_Analysis';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', ...
         'xScale', 'log', 'yScale', 'log');

% Plot FFT
plot(Sr(1:floor(n / 2)), rescale(PSD(1:floor(n / 2))), ...
     'color', graphColours(1), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1e-2; 1e3]);
ylim([1e-13; 1e2]);
tickData = [1e-1; 1e0; 1e1; 1e2];
xticks(tickData);
tickData = [1e-10; 1e-7; 1e-4; 1e-1];
yticks(tickData);
% xtickformat('%.1f');
% ytickformat('%.1f');
xlabel({'{$St_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{PSD}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');