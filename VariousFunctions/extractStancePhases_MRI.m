function [IC_L,IC_R,TO_L,TO_R] = extractStancePhases_MRI(path_GRF,threshold,side)

GRF = importdata(path_GRF);
GRF_L_vy = GRF.data(:,strcmp('     ground_force_vy',GRF.colheaders));
GRF_R_vy = GRF.data(:,strcmp('   1_ground_force_vy',GRF.colheaders));

switch side
    case 'LEFT'
        if isempty(GRF_L_vy)
            GRF_L_vy = GRF.data(:,strcmp('ground_force_vy',GRF.colheaders));
            GRF_R_vy = GRF.data(:,strcmp('1_ground_force_vy',GRF.colheaders));
        end
    case 'RIGHT'
        if isempty(GRF_R_vy)
            GRF_L_vy = GRF.data(:,strcmp('ground_force_vy',GRF.colheaders));
            GRF_R_vy = GRF.data(:,strcmp('1_ground_force_vy',GRF.colheaders));
        end
end

GRF_time = GRF.data(:,strcmp('time',GRF.colheaders));

if isempty(GRF_time)
    GRF_time = GRF.data(:,strcmp('                time',GRF.colheaders));
end
logic_L = GRF_L_vy >= threshold;
logic_R = GRF_R_vy >= threshold;

GRF_L_IC = find(diff(logic_L)==1) +1;
GRF_R_IC = find(diff(logic_R)==1) +1;

GRF_L_TO = (diff(logic_L)==-1);
GRF_R_TO = (diff(logic_R)==-1);

IC_L = GRF_time(GRF_L_IC);
IC_R = GRF_time(GRF_R_IC);

TO_L = GRF_time(GRF_L_TO);
TO_R = GRF_time(GRF_R_TO);

end


