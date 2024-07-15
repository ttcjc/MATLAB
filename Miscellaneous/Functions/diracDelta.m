%% Dirac Delta Function v1.0
% ----
% Generates a Unit Impulse, With a Value of Zero Everywhere Except at Zero 
% ----
% Usage: output = diracDelta(input, tol)
%
%        'input'   -> Any Numeric Value (Scalar, Vector or Array)
%        'tol'     -> Number of Decimal Places Used When Rounding to Provide a Tolerance for Zero Values


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function output = diracDelta(input, tol)
    
    output = (round(input, tol) == 0) * 1;

end