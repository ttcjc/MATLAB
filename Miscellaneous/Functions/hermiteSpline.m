%% Blah
% ----
% Blah
% ----
% Usage: Blah


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [x, y] = hermiteSpline(P1, P2, m1, m2, nP)

    % Set Interpolation Resolution
    t = linspace(0, 1, nP)';
    
    % Define Hermite Basis Functions
    H11 = (2 * t.^3) - (3 * t.^2) + 1;
    H12 = (- 2 * t.^3) + (3 * t.^2);
    H21 = t.^3 - (2 * t.^2) + t;
    H22 = t.^3 - t.^2;
    
    % Calculate Spline
    x = P1(1) + (P2(1) - P1(1)) * t;
    y = (H11 * P1(2)) + (H21 * m1 * (P2(1) - P1(1))) + (H12 * P2(2)) + (H22 * m2 * (P2(1) - P1(1)));
    
end