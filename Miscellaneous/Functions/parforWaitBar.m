function parforWaitBar(handle, iterations)
    
    persistent i wB N
    
    if nargin == 2
        % Initialise Wait Bar
        i = 0;
        wB = handle;
        N = iterations;
    else
        % Update Wait Bar
        if isvalid(wB)
            i = i + 1;
            waitbar((i / N), wB);
        end
        
    end
    
end
