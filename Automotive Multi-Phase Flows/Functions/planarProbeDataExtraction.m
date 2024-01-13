%% Probe Data Initialisation v2.1
% ----
% Extract Planar Data From Previously Collated Volume Probe Data
% ----
% Usage: planarData = planarProbeDataExtraction(probeData, spacePrecision, format);
%
%        'probeData'      -> Volumetric Velocity Data Collated Using 'initialiseProbeData'
%        'spacePrecision' -> Desired Rounding Precision
%        'format'         -> Data Extraction Format ['singlePlane' / 'multiPlane']


%% Changelog

% v1.0 - Initial Commit
    

%% Main Function

function planarData = planarProbeDataExtraction(probeData, spacePrecision, format)

    % Identify Desired Plane(s)
    if strcmp(format, 'singlePlane')
        volumeSlice = identifyVolumeSlices(probeData.positionGrid, spacePrecision, false);
    else
        volumeSlice = identifyVolumeSlices(probeData.positionGrid, spacePrecision, true);
    end
    
    planes = fieldnames(volumeSlice);
    
    for i = 1:height(planes)
        
        % Identify Data in Plane
        switch volumeSlice.(planes{i}).orientation

            case 'YZ'
                index = find((probeData.positionGrid(:,1) == volumeSlice.(planes{i}).xLims) & ...
                             ((probeData.positionGrid(:,2) >= volumeSlice.(planes{i}).yLims(1)) & ...
                              (probeData.positionGrid(:,2) <= volumeSlice.(planes{i}).yLims(2))) & ...
                             ((probeData.positionGrid(:,3) >= volumeSlice.(planes{i}).zLims(1)) & ...
                              (probeData.positionGrid(:,3) <= volumeSlice.(planes{i}).zLims(2))));
                
            case 'XZ'
                index = find(((probeData.positionGrid(:,1) >= volumeSlice.(planes{i}).xLims(1)) & ...
                              (probeData.positionGrid(:,1) <= volumeSlice.(planes{i}).xLims(2))) & ...
                             (probeData.positionGrid(:,2) == volumeSlice.(planes{i}).yLims) & ...
                             ((probeData.positionGrid(:,3) >= volumeSlice.(planes{i}).zLims(1)) & ...
                              (probeData.positionGrid(:,3) <= volumeSlice.(planes{i}).zLims(2))));
                
            case 'XY'
                index = find(((probeData.positionGrid(:,1) >= volumeSlice.(planes{i}).xLims(1)) & ...
                              (probeData.positionGrid(:,1) <= volumeSlice.(planes{i}).xLims(2))) & ...
                             ((probeData.positionGrid(:,2) >= volumeSlice.(planes{i}).yLims(1)) & ...
                              (probeData.positionGrid(:,2) <= volumeSlice.(planes{i}).yLims(2))) & ...
                             (probeData.positionGrid(:,3) == volumeSlice.(planes{i}).zLims));
                
        end 

        % Collate Planar Data
        planarData.(planes{i}).positionGrid = probeData.positionGrid(index,:);
        
        planarData.(planes{i}).time = probeData.time;
        
        planarData.(planes{i}).u.inst = cell(height(planarData.(planes{i}).time),1);
        planarData.(planes{i}).v.inst = planarData.(planes{i}).u.inst;
        planarData.(planes{i}).w.inst = planarData.(planes{i}).u.inst;
        
        for j = 1:height(planarData.(planes{i}).time)
            planarData.(planes{i}).u.inst{j} = probeData.u.inst{j}(index);
            planarData.(planes{i}).v.inst{j} = probeData.v.inst{j}(index);
            planarData.(planes{i}).w.inst{j} = probeData.w.inst{j}(index);
        end
        clear j;
        
    end
    clear i;
    
end