function sstruct = ForwardSimulation_SplineInputData(t,input)

numColPoints = length(t);
NMuscles = input.auxdata.NMuscles;

sstruct.Excitation = zeros(numColPoints,NMuscles);
sstruct.lMT = zeros(numColPoints,NMuscles);
sstruct.vMT = zeros(numColPoints,NMuscles);

for m = 1:NMuscles
     sstruct.EMG(:,m) = ppval(input.auxdata.EMGSpline(m),t);
     [sstruct.lMT(:,m),sstruct.vMT(:,m),~] = ...
         SplineEval_ppuval(input.auxdata.LMTSpline(m),t,1);
end
