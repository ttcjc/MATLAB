clc;
close all;


%% 

y1 = sin(linspace(0, tau, 1000)');
y2 = 2 * sin(linspace(0, tau, 1000)');
y3 = sin( 2 * (linspace(0, tau, 1000)'));
y4 = sin(linspace(0, tau, 1000)').^2;


%%

xLimsPlot = [0.3; 2.1];
yLimsPlot = [-0.3; 0.3];
zLimsPlot = [0; 0.5];
                
                
%%
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [1200, 500, 1590, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

plot(y1, 'lineWidth', 2);
plot(y2, 'lineWidth', 2);
plot(y3, 'lineWidth', 2);
plot(y4, 'lineWidth', 2);

title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
lightangle(0, 45);
axis on;
box on;
grid off;
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([zLimsPlot(1), zLimsPlot(2)]);
tickData = xLimsPlot(1):(diff(xLimsPlot) / 5):xLimsPlot(2);
xticks(tickData(2:5));
tickData = zLimsPlot(1):(diff(zLimsPlot) / 5):zLimsPlot(2);
yticks(tickData(2:5));
xtickformat('%+.2g');
ytickformat('%+.2g');
xlabel({'{$x_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$z_{_{\ell}}$}'}, 'interpreter', 'latex');
legend({'Uncoupled', ...
        'Coupled', ...
        'Reduced Injection Mass', ...
        'Increased Injection Arc'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', 'Test', '.png'], '-dpng', '-r300');

