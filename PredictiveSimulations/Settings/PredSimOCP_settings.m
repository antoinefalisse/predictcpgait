% settings describes the parameters used in the optimal control problem.
% settings(1): weight joint kinematics
% settings(2): weight ground reaction forces
% settings(3): weight ground reaction moments
% settings(4): weight joint kinetics
% settings(5): weight muscle activations
% settings(6): number of mesh intervals
% settings(7): NLP error tolerance: 1*10^(-settings(9)).
% settings(8): bounds (1: min-max, 2: envelope).
% settings(9): index trials
% settings(10): index muscles for which no FLV
% settings(11): number of contact spheres
% settings(12): number of synergies
% settings(13): weight on synergies
% settings(14): allowed deviation from synergy weights
% settings(15): allowed deviation from pelvis residuals
% settings(16): index MT-parameters
% settings(17): index initial guess
% settings(18): index periodicity
% settings(19): whether increasing Fmax
% settings(20): track dynamically consistent motion 
% settings(ww,21): weight metabolic energy
% settings(ww,22): exponent metabolic energy


settings = [               
    % Increase weight on accelerations + no synergies
    0,   2500 ,150,4,99,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %695 Best
    0,   2500 ,150,4,99,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %696
    0,   2500 ,150,4,99,1,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %697 Best
    0,   2500 ,150,4,99,1,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %698
    100, 2500 ,150,4,99,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %699 Best
    100, 2500 ,150,4,99,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %700
    % Add 4 synergies
    0,   2500 ,150,4,4,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %701 Best
    0,   2500 ,150,4,4,1,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %702 Best
    100, 2500 ,150,4,4,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %703    
    0,   2500 ,150,4,4,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %704
    0,   2500 ,150,4,4,1,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %705
    100, 2500 ,150,4,4,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %706    
    0,   2500 ,150,4,4,2,5,0.25,2,10,0,50,1,0.1,1000, 695,;     %707
    0,   2500 ,150,4,4,1,5,0.25,2,10,0,50,1,0.1,1000, 697,;     %708
    100, 2500 ,150,4,4,2,5,0.25,2,10,0,50,1,0.1,1000, 6;     %709 Best
    % Add 3 synergies
    0,   2500 ,150,4,3,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %710
    0,   2500 ,150,4,3,1,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %711 
    100, 2500 ,150,4,3,2,4,0.25,2,10,0,50,1,0.1,1000, 0,;       %712    
    0,   2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %713
    0,   2500 ,150,4,3,1,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %714
    100, 2500 ,150,4,3,2,3,0.25,2,10,0,50,1,0.1,1000, 0,;       %715 
    0,   2500 ,150,4,3,2,5,0.25,2,10,0,50,1,0.1,1000, 695,;     %716 Best
    0,   2500 ,150,4,3,1,5,0.25,2,10,0,50,1,0.1,1000, 697,;     %717
    100, 2500 ,150,4,3,2,5,0.25,2,10,0,50,1,0.1,1000, 6;     %718 Best
    0,   2500 ,150,4,3,2,5,0.25,2,10,0,50,1,0.1,1000, 701,;     %719
    0,   2500 ,150,4,3,1,5,0.25,2,10,0,50,1,0.1,1000, 702,;     %720 Best
    100, 2500 ,150,4,3,2,5,0.25,2,10,0,50,1,0.1,1000, 709,;     %721
    % Include spasticity
    0,   2500 ,150,4,99,2,4,0.25,2,10,1,50,1,0.1,1000, 0,;       %723
    0,   2500 ,150,4,99,2,3,0.25,2,10,1,50,1,0.1,1000, 0,;       %724
    0,   2500 ,150,4,99,2,5,0.25,2,10,1,50,1,0.1,1000, 695,;     %725 Best
    % Add 4 synergies to spasticity
    0,   2500 ,150,4,4,2,4,0.25,2,10,1,50,1,0.1,1000, 0,;       %726 Max (20,000)
    0,   2500 ,150,4,4,2,3,0.25,2,10,1,50,1,0.1,1000, 0,;       %727 
    0,   2500 ,150,4,4,2,5,0.25,2,10,1,50,1,0.1,1000, 701,;     %728 Max (lower than 729...)
    0,   2500 ,150,4,4,2,5,0.25,2,10,1,50,1,0.1,1000, 725,;     %729 Best
    
    
    
    
    
    
   
    ];