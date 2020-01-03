function sstruct = forceModel_SplineInputData(t,input,ms)

numColPoints    = length(t);
NMuscles_spas   = input.auxdata.NMuscles_spas;

sstruct.Ftilde = zeros(numColPoints,NMuscles_spas);
sstruct.dFtilde = zeros(numColPoints,NMuscles_spas);
sstruct.EMG     = zeros(numColPoints,NMuscles_spas);

for m = 1:NMuscles_spas
    sstruct.Ftilde(:,m)    = ppval(input.auxdata.FtildeSpline(ms).MS(m),t);
    sstruct.dFtilde(:,m)    = ppval(input.auxdata.dFtildeSpline(ms).MS(m),t);
    sstruct.EMG(:,m)        = ppval(input.auxdata.EMGSpline(ms).MS(m),t);
end
