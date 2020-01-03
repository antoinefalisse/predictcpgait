% This script generates a series of plots showing bounds(mpi).mp and initial guess(mpi).mp
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
suptitle_Fontsize = 16;
figure()
for i = 1:size(bounds(mpi).mp.QsQdots.lower,2)
    subplot(10,6,i)
    plot(1:N,bounds(mpi).mp.QsQdots.upper(:,i)',...
        'b--','linewidth',2);
    hold on
    plot(1:N,bounds(mpi).mp.QsQdots.lower(:,i)',...
        'r--','linewidth',2);
    plot(guess(mpi).mp.QsQdots(:,i),'k','linewidth',2);
end
s = suptitle('Qs and Qdots');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:size(bounds(mpi).mp.Qdotdots.lower,2)
    subplot(5,6,i)
    plot([1,N],[bounds(mpi).mp.Qdotdots.upper(:,i),bounds(mpi).mp.Qdotdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(mpi).mp.Qdotdots.lower(:,i),bounds(mpi).mp.Qdotdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(mpi).mp.Qdotdots(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of Qdots');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMact
    subplot(10,10,i)
    plot([1,N],[bounds(mpi).mp.a.upper(:,i),bounds(mpi).mp.a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(mpi).mp.a.lower(:,i),bounds(mpi).mp.a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(mpi).mp.a(:,i),'k','linewidth',2);
end
s = suptitle('Muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMact
    subplot(10,10,i)
    plot([1,N],[bounds(mpi).mp.vA.upper(:,i),bounds(mpi).mp.vA.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(mpi).mp.vA.lower(:,i),bounds(mpi).mp.vA.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(mpi).mp.vA(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscles_FLV
    subplot(10,10,i)
    plot([1,N],[bounds(mpi).mp.FTtilde.upper(:,i),bounds(mpi).mp.FTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(mpi).mp.FTtilde.lower(:,i),bounds(mpi).mp.FTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(mpi).mp.FTtilde(:,i),'k','linewidth',2);
end
s = suptitle('Muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscles_FLV
    subplot(10,10,i)
    plot([1,N],[bounds(mpi).mp.dFTtilde.upper(:,i),bounds(mpi).mp.dFTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(mpi).mp.dFTtilde.lower(:,i),bounds(mpi).mp.dFTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(mpi).mp.dFTtilde(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
% figure()
% for i = 1:nq.trunk
%     subplot(3,3,i)
%     plot([1,N],[bounds(mpi).mp.a_b.upper(:,i),bounds(mpi).mp.a_b.upper(:,i)],...
%         'b--','linewidth',2);
%     hold on
%     plot([1,N],[bounds(mpi).mp.a_b.lower(:,i),bounds(mpi).mp.a_b.lower(:,i)],...
%         'r--','linewidth',2);
%     plot(guess(mpi).mp.a_b(:,i),'k','linewidth',2);
% end
% s = suptitle('Trunk activations');
% set(s,'Fontsize',suptitle_Fontsize)
% figure()
% for i = 1:nq.trunk
%     subplot(3,3,i)
%     plot([1,N],[bounds(mpi).mp.e_b.upper(:,i),bounds(mpi).mp.e_b.upper(:,i)],...
%         'b--','linewidth',2);
%     hold on
%     plot([1,N],[bounds(mpi).mp.e_b.lower(:,i),bounds(mpi).mp.e_b.lower(:,i)],...
%         'r--','linewidth',2);
%     plot(guess(mpi).mp.e_b(:,i),'k','linewidth',2);
% end
% s = suptitle('Trunk excitations');
% set(s,'Fontsize',suptitle_Fontsize)

% %% Evalue passive forces
% QsQdots_nsc = guess(mpi).mp.QsQdots.*repmat(scaling(mpi).mp.QsQdots,N,1);
% for k = 1:N
% % Get muscle-tendon lengths, velocities, and moment arms
%         % Left leg
%         qin_l_g = [QsQdots_nsc(k,jointi.hip_flex.l*2-1),...
%             QsQdots_nsc(k,jointi.hip_add.l*2-1), ...
%             QsQdots_nsc(k,jointi.hip_rot.l*2-1), ...
%             QsQdots_nsc(k,jointi.knee.l*2-1), ...
%             QsQdots_nsc(k,jointi.ankle.l*2-1),...
%             QsQdots_nsc(k,jointi.subt.l*2-1)];  
%         qdotin_l_g = [QsQdots_nsc(k,jointi.hip_flex.l*2),...
%             QsQdots_nsc(k,jointi.hip_add.l*2),...
%             QsQdots_nsc(k,jointi.hip_rot.l*2),...
%             QsQdots_nsc(k,jointi.knee.l*2),...
%             QsQdots_nsc(k,jointi.ankle.l*2),...
%             QsQdots_nsc(k,jointi.subt.l*2)];  
%         [lMTk_l_g,vMTk_l_g,~] = f_lMT_vMT_dM_l(qin_l_g,qdotin_l_g);    
%         % Right leg
%         qin_r_g = [QsQdots_nsc(k,jointi.hip_flex.r*2-1),...
%             QsQdots_nsc(k,jointi.hip_add.r*2-1),...
%             QsQdots_nsc(k,jointi.hip_rot.r*2-1),...
%             QsQdots_nsc(k,jointi.knee.r*2-1),...
%             QsQdots_nsc(k,jointi.ankle.r*2-1),...
%             QsQdots_nsc(k,jointi.subt.r*2-1)];  
%         qdotin_r_g = [QsQdots_nsc(k,jointi.hip_flex.r*2),...
%             QsQdots_nsc(k,jointi.hip_add.r*2),...
%             QsQdots_nsc(k,jointi.hip_rot.r*2),...
%             QsQdots_nsc(k,jointi.knee.r*2),...
%             QsQdots_nsc(k,jointi.ankle.r*2),...
%             QsQdots_nsc(k,jointi.subt.r*2)];      
%         [lMTk_r_g,vMTk_r_g,~] = f_lMT_vMT_dM_r(qin_r_g,qdotin_r_g);
%         % Both legs
%         lMTk_lr_g = [full(lMTk_l_g);full(lMTk_r_g)];
%         vMTk_lr_g = [full(vMTk_l_g);full(vMTk_r_g)];   
%         % Get muscle-tendon forces and derive Hill-equilibrium
%         [~,FT,~,Fpass_t,~,~,~] = ...
%             f_forceEquilibrium_FtildeState_all(guess(mpi).mp.a(k,musi_FLV),...
%                 guess(mpi).mp.FTtilde(k,:).*scaling(mpi).mp.FTtilde,guess(mpi).mp.dFTtilde(k,:).*scaling(mpi).mp.dFTtilde,...
%                 lMTk_lr_g(musi_FLV),vMTk_lr_g(musi_FLV),tensions(musi_FLV)); 
%         Fpass_g(k,:) = full(Fpass_t)';  
%         FT_g(k,:) = full(FT)';     
% end

