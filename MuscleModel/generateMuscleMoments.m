function [M_muscles,M_knee,FT] = generateMuscleMoments(time,Parameters,act_side,lMT_side,MA_side,Fvparam,Faparam,ind_muscles,Kpe)

Parameters_muscles = Parameters(:,ind_muscles);
lMT_muscles = lMT_side(:,ind_muscles);
MA_muscles = MA_side(:,ind_muscles);
act_muscles = act_side(:,ind_muscles);

FT = zeros(size(lMT_muscles));

% Distinguish between vastii (VI,VL,VM), medial hamstrings (BFLH,SM,ST) and regarding Kpe
% Kpe(1) for muscles for which stiffness is maintained generic, indices [2 3 4 5 6 7]
% Kpe(2) for medial hamstrings, indices [1 8 9]
% Kpe(3) for vastii, indices [10 11 12]

for j = 1:length(ind_muscles)
    i = ind_muscles(j);
    if (1<i && i<8)
        Kpef = Kpe(1);
    elseif (9<i)
        Kpef = Kpe(3);
    else
        Kpef = Kpe(2);
    end
    pp_y = spline(time,lMT_muscles(:,j));
    [LMTg,vMTg,~] = SplineEval_ppuval(pp_y,time,1);
    [FT(:,j),~] = integrateContractionDynamics(time,act_muscles(:,j),LMTg,vMTg,Parameters_muscles(:,j),Fvparam,Faparam,Kpef,'implicit');
end

M_muscles = MA_muscles.*FT;
M_knee = sum(M_muscles,2);

end






