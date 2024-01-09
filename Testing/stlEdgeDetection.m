% https://uk.mathworks.com/matlabcentral/answers/467744-efficient-algorithm-for-plotting-edges-detect-in-a-triangular-mesh?s_tid=prof_contriblnk

run preamble;

[geometry, ~, ~, ~, ~, ~] = selectGeometry(geoLoc);

parts = fieldnames(geometry);

for i = 1:height(parts)
    
    % Calculate Delaunay Triangulation
    con = delaunay(geometry.(parts{i}).vertices(:,1), ...
                   geometry.(parts{i}).vertices(:,2), ...
                   geometry.(parts{i}).vertices(:,3));
    
    % Identify Mesh Edges
    M = [con(:,[1,2]); con(:,[2,3]); con(:,[3,1])];   % Find all edges in mesh, note internal edges are repeated
    [u, ~, n] = unique(sort(M,2), 'rows'); % determine uniqueness of edges
    counts = accumarray(n(:), 1);   % determine counts for each unique edge
    O = u((counts == 1),:); % extract edges that only occurred once
    I = O(:,[1,2])';
    [x0, y0, z0] = deal(geometry.(parts{i}).vertices(I,1), ...
                        geometry.(parts{i}).vertices(I,1), ...
                        geometry.(parts{i}).vertices(I,1));
    
    % Plot Results
    figure()
    hold on;
%     trisurf(con, geometry.(parts{i}).vertices(:,1), ...
%                  geometry.(parts{i}).vertices(:,2), ...
%                  geometry.(parts{i}).vertices(:,3))
    plot3(x0, y0, z0, 'color', 'r', 'lineWidth', 2)
    view([30,30])
end