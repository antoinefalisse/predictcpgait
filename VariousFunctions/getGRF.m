% This function returns the GRF given as input the path to the
% file containing the raw data

function GRF = getGRF(pathGRF)

GRFall = importdata(pathGRF);
GRF.time = GRFall.data(:,strcmp(GRFall.colheaders,{'time'}));
identifier = {'1_ground_force_v','ground_force_v'}; % Right leg: 1_ground_force_v; Left leg: ground_force_v
identifierp = {'1_ground_force_p','ground_force_p'}; % Right leg: 1_ground_force_p; Left leg: ground_force_p
identifierm = {'1_ground_torque_','ground_torque_'}; % Right leg: 1_ground_torque_; Left leg: ground_torque_
if isempty(GRF.time)
    GRF.time = GRFall.data(:,strcmp(GRFall.colheaders,{'                time'}));
    identifier = {'   1_ground_force_v','     ground_force_v'}; % Right leg: 1_ground_force_v; Left leg: ground_force_v
    identifierp = {'   1_ground_force_p','     ground_force_p'}; % Right leg: 1_ground_force_p; Left leg: ground_force_p
    identifierm = {'   1_ground_torque_','     ground_torque_'}; % Right leg: 1_ground_torque_; Left leg: ground_torque_
end
% Weird additional case for simulation computer
if sum(strcmp(GRFall.colheaders,[identifier{1},'x']))==0 && sum(strcmp(GRFall.colheaders,{'time'}))==1
    identifier = {'   1_ground_force_v','     ground_force_v'}; % Right leg: 1_ground_force_v; Left leg: ground_force_v
    identifierp = {'   1_ground_force_p','     ground_force_p'}; % Right leg: 1_ground_force_p; Left leg: ground_force_p
    identifierm = {'   1_ground_torque_','     ground_torque_'}; % Right leg: 1_ground_torque_; Left leg: ground_torque_
end
    
GRF.val.all(:,1) = GRF.time;
GRF.pos.all(:,1) = GRF.time;
GRF.Mcop.all(:,1) = GRF.time;
GRF.MorGF.all(:,1) = GRF.time;
axis = {'x','y','z'};
leg = {'r','l'};
count = 1;
for j = 1:length(leg)
    for i = 1:length(axis)    
        count = count + 1;
        GRF.val.(leg{j})(:,i) = GRFall.data(:,strcmp(GRFall.colheaders,[identifier{j},axis{i}]));
        GRF.val.all(:,count) = GRF.val.(leg{j})(:,i);
        GRF.pos.(leg{j})(:,i) = GRFall.data(:,strcmp(GRFall.colheaders,[identifierp{j},axis{i}]));
        GRF.pos.all(:,count) = GRF.pos.(leg{j})(:,i);
        GRF.Mcop.(leg{j})(:,i) = GRFall.data(:,strcmp(GRFall.colheaders,[identifierm{j},axis{i}]));
        GRF.Mcop.all(:,count) = GRF.Mcop.(leg{j})(:,i);
    end
end
GRF.Mcop.Y = GRF.Mcop.all(:,[3,6]); 

% Calculate moments wrt ground frame origin (following Gil's notes, also
% confirmed with OpenSim and Simbody functions in C++)
for j = 1:length(leg)
    GRF.MorGF.(leg{j})(:,1) =  GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,3) - GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,2);
    GRF.MorGF.(leg{j})(:,2) =  GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,1) - GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,3) + GRF.Mcop.(leg{j})(:,2);
    GRF.MorGF.(leg{j})(:,3) =  GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,2) - GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,1);
end
GRF.MorGF.all(:,2:4) = GRF.MorGF.r;
GRF.MorGF.all(:,5:7) = GRF.MorGF.l;
