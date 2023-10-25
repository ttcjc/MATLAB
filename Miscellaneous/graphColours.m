%% Graph Colours v1.0
% ----
% Loughborough University Brand Colours for Use in Graphs
% ----
% Usage: colour = colour = graphColours(n)
%        'n' -> Desired Colour (1-7)


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function colour = graphColours(n)

    colours = [
               '#4A1863'
               '#E5007D'
               '#22C4AC'
               '#FCC21D'
               '#F9762D'
               '#575756'
               '#1D1D1B'
              ];
          
    colour = colours(n,:);

end