%% Large Rainbow Colour Map v1.0
% ----
% A (Dreadful) Seven Point Color Spectrum From Blue to Cyan to Green to Yellow To Red to Purple to White
%  Featured in Tecplot
% ----
% Usage: map = largeRainbow(m)
%        'm' -> Desired Map Levels


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function map = largeRainbow(m)

    if nargin < 1
        f = get(groot,'CurrentFigure');
        
        if isempty(f)
            m = height(get(groot,'DefaultFigureColormap'));
        else
            m = height(f.Colormap);
        end
    
    end
    
    values = [
              0 0 255
              0 255 255
              0 255 0
              255 255 0
              255 0 0
              255 0 255
              255 255 255
             ] / 255;
    
    P = height(values);
    map = interp1(1:P, values, linspace(1,P,m), 'linear');

end