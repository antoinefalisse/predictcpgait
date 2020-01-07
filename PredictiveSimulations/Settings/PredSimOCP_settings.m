% settings contains parameters used in the optimal control problem.

% settings(ww,1);	% weight TD kinematics tracking
% settings(ww,2);   % weight muscle activations (fatigue)
% settings(ww,3);   % number of mesh intervals
% settings(ww,4);   % NLP error tolerance: 1*10^(-tol_ipopt)
% settings(ww,5);   % number of synergies per leg
% settings(ww,6);   % MT-parameters (generic-1 or personalized-2)
% settings(ww,7);   % initial guess identifier: 1-CP walking pattern, 2-TD
%   walking pattern, 3-movement simulated using the row in settings given by
%   settings(ww,16) (e.g., case 15 uses case 5 as initial guess)
% settings(ww,8);   % weight metabolic energy rate
% settings(ww,9);   % exponent metabolic energy rate
% settings(ww,10);  % exponent muscle activations
% settings(ww,11);  % whether to take spasticity into account (yes-1, no-0)
% settings(ww,12);  % weight joint accelerations
% settings(ww,13);  % imposed gait speed (m s-1)
% settings(ww,14);  % weight trunk excitations
% settings(ww,15);  % weight passive joint torques
% settings(ww,16);  % initial guess case (see comment for settings(ww,7))

settings = [     
    %% no synergies %%
    % Personalized MT-parameters (guesses 2 & 1)
    0,   2500 ,150,4,99,2,2,0.25,2,10,0,50,1,0.1,1000, 0;	%1 Best optimal cost among guesses
    0,   2500 ,150,4,99,2,1,0.25,2,10,0,50,1,0.1,1000, 0;	%2
    % Generic MT-parameters (guesses 2 & 1)
    0,   2500 ,150,4,99,1,2,0.25,2,10,0,50,1,0.1,1000, 0;	%3 Best optimal cost among guesses
    0,   2500 ,150,4,99,1,1,0.25,2,10,0,50,1,0.1,1000, 0;	%4
    % Personalized MT-parameters + tracking of TD kinematics (guesses 2 & 1)
    100, 2500 ,150,4,99,2,2,0.25,2,10,0,50,1,0.1,1000, 0;  	%5 Best optimal cost among guesses
    100, 2500 ,150,4,99,2,1,0.25,2,10,0,50,1,0.1,1000, 0;  	%6
    %% 4 synergies %%
    % Personalized MT-parameters (guess 2)
    0,   2500 ,150,4,4,2,2,0.25,2,10,0,50,1,0.1,1000, 0;  	%7 Best optimal cost among guesses
    % Generic MT-parameters (guess 2)
    0,   2500 ,150,4,4,1,2,0.25,2,10,0,50,1,0.1,1000, 0;  	%8 Best optimal cost among guesses
    % Personalized MT-parameters + tracking of TD kinematics (guess 2)
    100, 2500 ,150,4,4,2,2,0.25,2,10,0,50,1,0.1,1000, 0;   	%9    
    % Personalized MT-parameters (guess 1)
    0,   2500 ,150,4,4,2,1,0.25,2,10,0,50,1,0.1,1000, 0;  	%10
    % Generic MT-parameters (guess 1)
    0,   2500 ,150,4,4,1,1,0.25,2,10,0,50,1,0.1,1000, 0;   	%11
    % Personalized MT-parameters + tracking of TD kinematics (guess 1)
    100, 2500 ,150,4,4,2,1,0.25,2,10,0,50,1,0.1,1000, 0; 	%12 
    % Personalized MT-parameters (guess 3)
    0,   2500 ,150,4,4,2,3,0.25,2,10,0,50,1,0.1,1000, 1; 	%13
    % Generic MT-parameters (guess 3)
    0,   2500 ,150,4,4,1,3,0.25,2,10,0,50,1,0.1,1000, 3; 	%14
    % Personalized MT-parameters + tracking of TD kinematics (guess 3)
    100, 2500 ,150,4,4,2,3,0.25,2,10,0,50,1,0.1,1000, 5;  	%15 Best optimal cost among guesses
    %% 3 synergies %%
    % Personalized MT-parameters (guess 2)
    0,   2500 ,150,4,3,2,2,0.25,2,10,0,50,1,0.1,1000, 0,;	%16
    % Generic MT-parameters (guess 2)
    0,   2500 ,150,4,3,1,2,0.25,2,10,0,50,1,0.1,1000, 0,;	%17 
    % Personalized MT-parameters + tracking of TD kinematics (guess 2)
    100, 2500 ,150,4,3,2,2,0.25,2,10,0,50,1,0.1,1000, 0,;  	%18    
    % Personalized MT-parameters (guess 1)
    0,   2500 ,150,4,3,2,1,0.25,2,10,0,50,1,0.1,1000, 0,;  	%19
    % Generic MT-parameters (guess 1)
    0,   2500 ,150,4,3,1,1,0.25,2,10,0,50,1,0.1,1000, 0,; 	%20
    % Personalized MT-parameters + tracking of TD kinematics (guess 1)
    100, 2500 ,150,4,3,2,1,0.25,2,10,0,50,1,0.1,1000, 0,;	%21 
    % Personalized MT-parameters (guess 3a)
    0,   2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 1; 	%22 Best optimal cost among guesses
    % Generic MT-parameters (guess 3a)
    0,   2500 ,150,4,3,1,3,0.25,2,10,0,50,1,0.1,1000, 3;  	%23
    % Personalized MT-parameters + tracking of TD kinematics (guess 3a)
    100, 2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 5;  	%24 Best optimal cost among guesses
    % Personalized MT-parameters (guess 3b)
    0,   2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 7;   	%25
    % Generic MT-parameters (guess 3b)
    0,   2500 ,150,4,3,1,3,0.25,2,10,0,50,1,0.1,1000, 8;  	%26 Best optimal cost among guesses
    % Personalized MT-parameters + tracking of TD kinematics (guess 3b)
    100, 2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 15; 	%27
    %% spasticity but no synergies %% 
    % Personalized MT-parameters (guess 2)
    0,   2500 ,150,4,99,2,2,0.25,2,10,1,50,1,0.1,1000, 0; 	%28
    % Personalized MT-parameters (guess 1)
    0,   2500 ,150,4,99,2,1,0.25,2,10,1,50,1,0.1,1000, 0; 	%29
    % Personalized MT-parameters (guess 3)
    0,   2500 ,150,4,99,2,3,0.25,2,10,1,50,1,0.1,1000, 1; 	%30 Best optimal cost among guesses
    %% spasticity + 4 synergies %%
    % Personalized MT-parameters (guess 2)
    0,   2500 ,150,4,4,2,2,0.25,2,10,1,50,1,0.1,1000, 0; 	%31
    % Personalized MT-parameters (guess 1)
    0,   2500 ,150,4,4,2,1,0.25,2,10,1,50,1,0.1,1000, 0;  	%32 
    % Personalized MT-parameters (guess 3a)
    0,   2500 ,150,4,4,2,3,0.25,2,10,1,50,1,0.1,1000, 7;   	%33
    % Personalized MT-parameters (guess 3b)
    0,   2500 ,150,4,4,2,3,0.25,2,10,1,50,1,0.1,1000, 30;  	%34 Best optimal cost among guesses
    ];
