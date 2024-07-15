%% Two-Dimensional Smoothing Function
% ----
% Smooths Two-Dimensional Data by Calculating the Spatial Average Within User-Defined Rectangles
% Centred on Each Data Point. 'NaN' Values Are Preserved and Do Not Impact the Average
% ----
% Usage: outputMatrix = smooth2(inputMatrix, windowSize)
%
%        'inputMatrix' -> Numeric Data in the Form of a Two-Dimensional 'ndgrid'
%        'windowSize'  -> Scalar or Two-Element Vector


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function outputMatrix = smooth2(inputMatrix, windowSize)
    
    % Confirm Validity of 'inputMatrix'
    if ~isnumeric(inputMatrix) || (numel(size(inputMatrix)) ~= 2)
        error('''inputMatrix'' Must Be Numeric Data in the Form of a Two-Dimensional ''ndgrid''');
    end
    
    % Confirm Validity of 'windowSize'
    if ~isnumeric(windowSize) || numel(windowSize) > 2
        error('Window Size Must Be a Scalar or a Two-Element Vector');
    end
    
    % Temporarily Half Window Size To Ensure Compatibility With 'spdiags'
    if numel(windowSize) == 2
        windowR = floor(windowSize(1) / 2);
        windowC = floor(windowSize(2) / 2);
    else
        windowR = floor(windowSize / 2);
        windowC = floor(windowSize / 2);
    end
    
    % Generate Smoothing Matrices
    [R, C] = size(inputMatrix);
    smoothR = spdiags(ones(R, ((2 * windowR) + 1)), (-windowR:windowR), R, R);
    smoothC = spdiags(ones(C, ((2 * windowC) + 1)), (-windowC:windowC), C, C);
    
    % Temporarily Remove NaN Values From 'inputMatrix'
    index = isnan(inputMatrix);
    inputMatrix(index) = 0;
    
    % Calculate Number of Valid (Numeric) Values Present in Each Window
    nValid = smoothR * (~index) * smoothC;
    nValid(index) = NaN;
    
    % Calculate Sum of Values in Each Window
    sumOfWindows = smoothR * inputMatrix * smoothC;
    
    % Calculate Spatial Average in Each Window
    outputMatrix = sumOfWindows ./ nValid;
    
end