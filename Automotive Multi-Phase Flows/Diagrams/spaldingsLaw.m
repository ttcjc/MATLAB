run preamble;


%%

channelFlow = [
               0.39775099, 0.41853505
               1.5927147, 1.6150399
               3.5804355, 3.5858682
               6.3777093, 6.0150847
               9.8887390, 8.4609338
               14.209386, 10.463720
               19.350229, 11.935120
               25.424230, 12.981467
               39.244642, 14.365335
               47.572575, 14.844690
               56.391646, 15.235378
               66.249959, 15.590598
               76.450781, 15.945681
               87.045701, 16.265263
               99.553462, 16.584879
               111.33898, 16.904323
               123.96368, 17.152937
               138.01989, 17.401550
               152.98374, 17.685527
               167.30849, 17.934002
               182.15749, 18.164745
               198.32437, 18.413186
               214.96191, 18.643893
               231.95477, 18.874567
               249.17327, 19.087507
               267.66995, 19.335845
               284.97745, 19.513317
               304.76495, 19.708524
               323.02212, 19.903661
               343.90869, 20.098833
               362.88312, 20.258538
               386.34709, 20.436011
               407.66297, 20.578016
               428.23410, 20.737687
               453.88778, 20.897426
               483.23611, 21.057200
               523.77504, 21.199412
               572.81910, 21.323994
              ];


%%

yPlusLinear = (0.01:0.01:15)';
uPlusLinear = yPlusLinear;

yPlusLog = (5:0.01:2000)';
uPlusLog = (1 / 0.41) * log(9 * yPlusLog);

uPlusSpalding = (0.01:0.01:25)';
yPlusSpalding = uPlusSpalding + 0.1108 * (exp(0.4 * uPlusSpalding) - 1 - (0.4 * uPlusSpalding) - ((0.4 * uPlusSpalding).^2 / factorial(2)) - ((0.4 * uPlusSpalding).^3 / factorial(3)) - ((0.4 * uPlusSpalding).^4 / factorial(4)));


%% 

% Initialise Figure
fig = fig + 1;
figName = 'Spaldings_Law';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', 'xScale', 'log');

% Plot Profiles
plot(channelFlow(:,1), channelFlow(:,2), 'color', graphColours(7), 'marker', 'o', ...
                                                                   'markerSize', 5, ...
                                                                   'lineStyle', 'none', ...
                                                                   'lineWidth', 2);
plot(yPlusSpalding, uPlusSpalding, 'color', graphColours(1), 'lineWidth', 2);
plot(yPlusLinear, uPlusLinear, 'color', graphColours(2),'lineStyle', ':', 'lineWidth', 2);
plot(yPlusLog, uPlusLog, 'color', graphColours(3),'lineStyle', ':', 'lineWidth', 2);

% Plot Regions
lineHandle = xline(5, 'alpha', 1, ...
                      'lineStyle', ':', ...
                      'lineWidth', 2, ...
                      'label', 'Laminar Sublayer $\quad$', ...
                      'labelHorizontalAlignment', 'Left', ...
                      'labelVerticalAlignment', 'Top');
lineHandle.Interpreter = 'latex';
lineHandle.FontSize = 18;
clear lineHandle;

lineHandle = xline(30, 'alpha', 1, ...
                       'lineStyle', ':', ...
                       'lineWidth', 2, ...
                       'label', '$\quad$ Buffer Region', ...
                       'labelHorizontalAlignment', 'Left', ...
                       'labelVerticalAlignment', 'Bottom');
lineHandle.Interpreter = 'latex';
lineHandle.FontSize = 18;
clear lineHandle;

lineHandle = xline(350, 'alpha', 1, ...
                        'lineStyle', ':', ...
                        'lineWidth', 2, ...
                        'label', '$\quad$ Inertial Sublayer', ...
                        'labelHorizontalAlignment', 'Left', ...
                        'labelVerticalAlignment', 'Bottom');
lineHandle.Interpreter = 'latex';
lineHandle.FontSize = 18;
clear lineHandle;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1e-1; 1e3]);
ylim([0; 25]);
tickData = [1e0; 1e1; 1e2];
xticks(tickData);
tickData = (5:5:20);
yticks(tickData);
xlabel({'{$y^{+}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$u^{+}$}'}, 'interpreter', 'latex');
legend({'DNS', ...
        'Spalding''s Law', ...
        '$u^{+}\,=\,y^{+}$', ...
        '$u^{+}\,=\,\frac{1}{\kappa}\,\ln(E\,y^{+})$'}, 'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
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