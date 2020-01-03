% This script returns the indices of the muscles actuating the joints, for
% use with the moment arms as well as the indices of the muscles in the
% vector containing all muscles

% order of dof: hip flex, hip add, hip rot, knee flex, ankle flex, sub

function [Indmusi,mai] = getMomentArmIndices(muscleNames,muscle_spanning_joint_INFO)

NMuscle = length(muscleNames)*2;

for i = 1:length(muscleNames)
    Indmusi.(muscleNames{i}(1:end-2)).l = find(strcmp(muscleNames,muscleNames{i}));
    Indmusi.(muscleNames{i}(1:end-2)).r = Indmusi.(muscleNames{i}(1:end-2)).l + NMuscle/2;
    
    for j = 1:size(muscle_spanning_joint_INFO,2)
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus.l(1,i) = Indmusi.(muscleNames{i}(1:end-2)).l;
           mai(j).mus.r(1,i) = Indmusi.(muscleNames{i}(1:end-2)).r;
        end        
    end    
end

for j = 1:size(muscle_spanning_joint_INFO,2)
    mai(j).mus.l(mai(j).mus.l == 0) = [];
    mai(j).mus.r(mai(j).mus.r == 0) = [];
end
end
