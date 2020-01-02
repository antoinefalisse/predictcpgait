function phaseout = Wrap4For_ForwardSimulation_cont(input)

global splinestruct

if isempty(splinestruct) || size(splinestruct.EMG,1) ~= length(input.phase.time) 
    splinestruct = ForwardSimulation_SplineInputData(input.phase.time,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = ForwardSimulation_cont(input);