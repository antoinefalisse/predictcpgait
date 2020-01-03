% Get muscle names/indices for full model/knee/knee extensors

Mall = {'glut_max1_r','glut_max2_r','glut_max3_r','glut_med1_r',...
    'glut_med2_r','glut_med3_r','glut_min1_r','glut_min2_r',...
    'glut_min3_r','add_long_r','add_brev_r','add_mag1_r','add_mag2_r',...
    'add_mag3_r','pectineus_r','iliacus_r','psoas_r','quad_fem_r',...
    'gemellus_r','piri_r','TFL_r','gracilis_r','semimem_r','semiten_r',...
    'bi_fem_lh_r','bi_fem_sh_r','sartorius_r','rectus_fem_r',...
    'vas_med_r','vas_int_r','vas_lat_r','gas_med_r',...
    'gas_lat_r','soleus_r','tib_post_r','tib_ant_r','ext_dig_r',...
    'ext_hal_r','flex_dig_r','flex_hal_r','per_brev_r','per_long_r',...
    'per_tert_r',...
    'glut_max1_l','glut_max2_l','glut_max3_l','glut_med1_l',...
    'glut_med2_l','glut_med3_l','glut_min1_l','glut_min2_l',...
    'glut_min3_l','add_long_l','add_brev_l','add_mag1_l','add_mag2_l',...
    'add_mag3_l','pectineus_l','iliacus_l','psoas_l','quad_fem_l',...
    'gemellus_l','piri_l','TFL_l','gracilis_l','semimem_l','semiten_l',...
    'bi_fem_lh_l','bi_fem_sh_l','sartorius_l','rectus_fem_l',...
    'vas_med_l','vas_int_l','vas_lat_l','gas_med_l',...
    'gas_lat_l','soleus_l','tib_post_l','tib_ant_l','ext_dig_l',...
    'ext_hal_l','flex_dig_l','flex_hal_l','per_brev_l','per_long_l',...
    'per_tert_l'};

% Knee muscles (2*12)
MKnee = {'bi_fem_lh_r'; 'bi_fem_sh_r'; 'gracilis_r'; 'gas_lat_r';...
    'gas_med_r'; 'rectus_fem_r'; 'sartorius_r'; 'semimem_r'; 'semiten_r';...
    'vas_int_r'; 'vas_lat_r'; 'vas_med_r'; 'bi_fem_lh_l';'bi_fem_sh_l';...
    'gracilis_l'; 'gas_lat_l'; 'gas_med_l'; 'rectus_fem_l'; 'sartorius_l';...
    'semimem_l'; 'semiten_l'; 'vas_int_l'; 'vas_lat_l'; 'vas_med_l'};

MAnkle = {'gas_med_r';'gas_lat_r';'soleus_r';'tib_post_r';'flex_dig_r';...
    'flex_hal_r';'tib_ant_r';'per_brev_r';'per_long_r';'per_tert_r';...
    'ext_dig_r';'ext_hal_r';'gas_med_l';'gas_lat_l';'soleus_l';...
    'tib_post_l';'flex_dig_l';'flex_hal_l';'tib_ant_l';'per_brev_l';...
    'per_long_l';'per_tert_l';'ext_dig_l';'ext_hal_l'};

% Knee extensors (2*4)
MKextensors = {'rectus_fem_r';'vas_int_r'; 'vas_lat_r'; 'vas_med_r';...
    'rectus_fem_l';'vas_int_l'; 'vas_lat_l'; 'vas_med_l'};

MAextensors = {'ext_dig_r';'ext_hal_r';'per_tert_r';'tib_ant_r';...
    'ext_dig_l';'ext_hal_l';'per_tert_l';'tib_ant_l'};

% Hamstrings (2*3)
MHamstrings = {'bi_fem_lh_r','semimem_r','semiten_r','bi_fem_lh_l',...
    'semimem_l','semiten_l'};

% Vastii (2*3)
MVastii = {'vas_int_r', 'vas_lat_r', 'vas_med_r',...
    'vas_int_l','vas_lat_l','vas_med_l'};

% Rectus femoris (2*1)
MRectus = {'rectus_fem_r','rectus_fem_l'};

% Number muscles per side
NMall_side = length(Mall)/2;
NMKnee_side = length(MKnee)/2;              
NMAnkle_side = length(MAnkle)/2;
NMKextensors_side = length(MKextensors)/2;  
NMAextensors_side = length(MAextensors)/2;
NMHamstrings_side = length(MHamstrings)/2;
NMVastii_side = length(MVastii)/2;
NMRectus_side = length(MRectus)/2;

% Select muscles as a function of side
switch allsegments(1).side_name
    case 'RIGHT'
        letter_side = 'r';
        vec_Mall_side = 1:NMall_side;
        vec_MKnee_side = 1:NMKnee_side;                 
        vec_MAnkle_side = 1:NMAnkle_side;                 
        vec_MKkextensors_side = 1:NMKextensors_side;    
        vec_MAextensors_side = 1:NMAextensors_side;    
        vec_MHamstrings_side = 1:NMHamstrings_side;
        vec_MVastii_side = 1:NMVastii_side;
        vec_MRectus_side = 1:NMRectus_side;

    case 'LEFT'
        letter_side = 'l';
        vec_Mall_side = NMall_side+1:2*NMall_side;
        vec_MKnee_side = NMKnee_side+1:2*NMKnee_side;                       
        vec_MAnkle_side = NMAnkle_side+1:2*NMAnkle_side;
        vec_MKkextensors_side = NMKextensors_side+1:2*NMKextensors_side;    
        vec_MAextensors_side = NMAextensors_side+1:2*NMAextensors_side;        
        vec_MHamstrings_side = NMHamstrings_side+1:2*NMHamstrings_side;        
        vec_MVastii_side = NMVastii_side+1:2*NMVastii_side;  
        vec_MRectus_side = NMRectus_side+1:2*NMRectus_side;  
end

Mall_side = Mall(vec_Mall_side);
MKnee_side = MKnee(vec_MKnee_side);                     
MAnkle_side = MAnkle(vec_MAnkle_side);                             
MKextensors_side = MKextensors(vec_MKkextensors_side);  
MAextensors_side = MAextensors(vec_MAextensors_side);
MHamstrings_side = MHamstrings(vec_MHamstrings_side);
MVastii_side = MVastii(vec_MVastii_side);
MRectus_side = MRectus(vec_MRectus_side);

% Indices MKnee_side in Mall_side
ind_MKnee = zeros(1,length(MKnee_side));
for i = 1:length(MKnee_side)
    ind_MKnee(i) = find(strcmp(Mall_side,MKnee_side{i})==1);
end

% Indices MAnkle_side in Mall_side
ind_MAnkle = zeros(1,length(MAnkle_side));
for i = 1:length(MAnkle_side)
    ind_MAnkle(i) = find(strcmp(Mall_side,MAnkle_side{i})==1);
end

% Indices kMKextensors_side in Mall_side
ind_MKextensors = zeros(1,length(MKextensors_side));
for i = 1:length(MKextensors_side)
    ind_MKextensors(i) = find(strcmp(Mall_side,MKextensors_side{i})==1);
end

% Indices kMAextensors_side in Mall_side
ind_MAextensors = zeros(1,length(MAextensors_side));
for i = 1:length(MAextensors_side)
    ind_MAextensors(i) = find(strcmp(Mall_side,MAextensors_side{i})==1);
end

% Indices MHamstrings in MMall_side
ind_MHamstrings = zeros(1,length(MHamstrings_side));
for i = 1:length(MHamstrings_side)
    ind_MHamstrings(i) = find(strcmp(Mall_side,MHamstrings_side{i})==1);
end

% Indices MVastii in MMall_side
ind_MVastii = zeros(1,length(MVastii_side));
for i = 1:length(MVastii_side)
    ind_MVastii(i) = find(strcmp(Mall_side,MVastii_side{i})==1);
end

% Indices MRectus in MMall_side
ind_MRectus = zeros(1,length(MRectus_side));
for i = 1:length(MRectus_side)
    ind_MRectus(i) = find(strcmp(Mall_side,MRectus_side{i})==1);
end

% Indices MHamstrings in MKnee_side
ind_MHamstringsKnee = zeros(1,length(MHamstrings_side));
for i = 1:length(MHamstrings_side)
    ind_MHamstringsKnee(i)=find(strcmp(MKnee_side,MHamstrings_side{i})==1);
end

% Indices MVastii in MKnee_side
ind_MVastiiKnee = zeros(1,length(MVastii_side));
for i = 1:length(MVastii_side)
    ind_MVastiiKnee(i) = find(strcmp(MKnee_side,MVastii_side{i})==1);
end
