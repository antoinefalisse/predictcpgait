function phaseout = forceModel_contGrdWrap(input)

global splinestruct

Nsegment = input.auxdata.Nsegment;

for ms = 1:Nsegment

    if isempty(splinestruct(ms).MS) || size(splinestruct(ms).MS.EMG,1) ~= length(input.phase(ms).time.f)    
        splinestruct(ms).MS = forceModel_SplineInputData(input.phase(ms).time.f,input,ms);
    end

    input.auxdata.splinestruct(ms).MS = splinestruct(ms).MS;    
end

phaseout = forceModel_contADiGatorGrd(input);