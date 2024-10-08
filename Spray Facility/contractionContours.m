run preamble;

%%%%

nP = 1e3;

Din = 2;

Dout = 1.15;

L = Din;

w1 = (Din / 2) / L;

w2 = (Dout / 2) / L;

X = (0.41 * L) / L;

x = linspace(0, 1, nP)';

%%%%

n1 = 4;

yUp = w1 - ((w1 - w2) * ((x.^(n1)) / (X^(n1 - 1))));

yUp_grad = gradient(yUp, mean(diff(x)));
 
%%%%

n2 = 6;

yDown = w2 - ((w2 - w1) * (((1 - x).^(n2)) / ((1 - X)^(n2 - 1))));

yDown_grad = gradient(yDown, mean(diff(x)));

%%%%

X_index = find((x <= X), 1, 'last');

[~, yUp_index] = min(yUp_grad(1:(X_index - 100)));

[~, yDown_index] = min(abs(yDown_grad - yUp_grad(yUp_index)));

%%%%

[xTransition, yTransition] = hermiteSpline([x(yUp_index), yUp(yUp_index)], ...
                                           [x(yDown_index), yDown(yDown_index)], ...
                                           yUp_grad(yUp_index), yDown_grad(yDown_index), ...
                                           round((x(yDown_index) - x(yUp_index)) * nP));

contractionProfile = [[x(1:yUp_index); xTransition; x(yDown_index:end)], ...
                      [yUp(1:yUp_index); yTransition; yDown(yDown_index:end)]] .* L;
contractionProfile = [contractionProfile, zeros(height(contractionProfile),1)];
contractionProfile = unique(contractionProfile, 'rows');

%%%%

% Initialise Figure
fig = fig + 1;
figName = 'Contour Gradients';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Contour Gradients
plot((x .* L), yUp_grad, 'color', graphColours(1), 'lineWidth', 2);
plot((x .* L), yDown_grad, 'color', graphColours(2), 'lineWidth', 2);
xline((X * L), 'alpha', 1, 'lineStyle', '--', 'lineWidth', 2);
plot((x(yUp_index) * L), yUp_grad(yUp_index), 'color', graphColours(1), 'lineWidth', 2, ...
                                              'lineStyle', 'none', 'marker', 'o', 'markerSize', 10);
plot((x(yDown_index) * L), yDown_grad(yDown_index), 'color', graphColours(2), 'lineWidth', 2, ...
                                                    'lineStyle', 'none', 'marker', 'o', 'markerSize', 10);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
axis padded
tickData = [];
xticks(tickData);
tickData = [];
yticks(tickData);
xlabel({'{$x$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$dy\,/\,dx$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

%%%%

% Initialise Figure
fig = fig + 1;
figName = 'Contraction Profile';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Contraction Profile
plot((x(1:yUp_index) .* L), (yUp(1:yUp_index) .* L), 'color', graphColours(1), 'lineWidth', 2);
plot((x(yDown_index:end) .* L), (yDown(yDown_index:end) .* L), 'color', graphColours(2), 'lineWidth', 2);
plot((xTransition .* L), (yTransition .* L), 'color', graphColours(3), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
axis padded
tickData = [];
xticks(tickData);
tickData = [];
yticks(tickData);
xlabel({'{$x$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$y$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

%%%%

contractionProfile = contractionProfile * 1000;
contractionProfile(:,1) = contractionProfile(:,1) - max(contractionProfile(:,1));
contractionProfile(:,3) = 575;
% writematrix(contractionProfile, '/home/lunet/ttcjc/Mount/OneDrive/04 - RA/Spray Facility/Generic Tunnel Design/contractionProfile.dat');