% Get muscle-tendon lengths and moment arms

function [lMT,MA] = getlMTMA(muscleNames,lMT_raw,MA_raw)

lMT = zeros(size(lMT_raw.data,1),length(muscleNames));
MA = zeros(size(MA_raw.data,1),length(muscleNames));

for i = 1:length(muscleNames)
    lMT(:,i) = lMT_raw.data(:,strcmp(lMT_raw.colheaders,muscleNames(i)));
    MA(:,i) = MA_raw.data(:,strcmp(MA_raw.colheaders,muscleNames(i)));
end

end
