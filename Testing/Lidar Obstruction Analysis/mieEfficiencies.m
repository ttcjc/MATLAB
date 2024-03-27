run preamble;


%% Initialise

n_a = 1.00027451; % Refractive Index of Air
n_w = complex(1.32352, 5.15e-7); % Refractive Index of Water

c = 299792458 / n_a; % Speed of Light in Air

k = (tau / 905e-9) * n_a; % Wave Number in Air


%% Calculate Mie Efficiencies

D = logspace(-8, -2, 10000)';

Qe_range = zeros([height(D),1]);
Qb_range = zeros([height(D),1]);

for i = 1:height(D)
    x = k * (D(i) / 2); % Size Parameter
    
    mieOut = mieScattering(n_w, x);
    
    Qe_range(i) = mieOut(1);
    Qb_range(i) = mieOut(4);
end
clear i mieOut;


%% Plot Mie Efficiencies

% Initialise Figure
fig = 1;
figName = 'Mie_Efficiencies';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top', 'xScale', 'log', 'yScale', 'log');

% Plot Efficiencies
plot((D * 1e6), Qb_range, 'color', graphColours(1), 'lineWidth', 2);
plot((D * 1e6), Qe_range, 'color', graphColours(2), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([1e-3; 1e5]);
ylim([1e-8; 1e2]);
tickData = [1e-1; 1e1; 1e3];
xticks(tickData);
tickData = [1e-6; 1e-4; 1e-2; 1e0];
yticks(tickData);
xlabel({'{$D_{_{p}}$ $(\mu m)$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Efficiency}'}, 'interpreter', 'latex');
legend({'{$Q_{_{b}}$}', '{$Q_{_{e}}$}'}, 'location', 'northWest', ...
                                 'orientation', 'vertical', ...
                                 'interpreter', 'latex', ...
                                 'fontSize', 18, ...
                                 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');