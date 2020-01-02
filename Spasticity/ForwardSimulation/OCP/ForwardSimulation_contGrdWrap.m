function phaseout = ForwardSimulation_contGrdWrap(input)

global splinestruct

if isempty(splinestruct) || size(splinestruct.EMG,1) ~= length(input.phase.time.f)
    splinestruct = ForwardSimulation_SplineInputData(input.phase.time.f,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = ForwardSimulation_contADiGatorGrd(input);