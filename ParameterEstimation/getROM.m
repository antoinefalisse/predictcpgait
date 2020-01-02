function [LMTout,PasMom,AllMuscles]=getROM(modelName,Values,Muscles,LUT,USEPASSIVEMOMENTS,Side)
% 
% computes the range of muscle lengths according to the data from the clinical exam 
% modelName..............string with the path and name of the osim model
% Values.................values of the range of motion and passive moments coming from the GUI
% Muscles................names of the muscles to be measured
% LUT....................look-up table to pass from the names in the GUI to the DoF in the model (optional -> [])
%                        {{~Name on GUI~       ~name in model~       ~name additional dof         value additional dof~ ...}}
% USEPASSIVEMOMENTS......whether or not to use information about the passive moments to tune muscle parameters (true/false)

% define the names and the default positions of the joints for each of the measurements
if isempty(LUT)
    LUT.H   =   {{'Flexion',            'hip_flex',     'knee_flex',    90}
                 {'Extension',          'hip_flex',     'knee_flex',    0}
                 {'Abduction 0',        'hip_add',      'knee_flex',    0,      'hip_flex',     0}
                 {'Abduction 90',       'hip_add',      'knee_flex',    90,     'hip_flex',     45}
                 {'Adduction',          'hip_add'}
                 {'Int Rot Sup',        'hip_rot',      'hip_flex',     90,     'knee_flex',    90}
                 {'Int Rot Pro',        'hip_rot',      'hip_flex',     0,      'knee_flex',    90}
                 {'Ext Rot Sup',        'hip_rot',      'hip_flex',     90,     'knee_flex',    90}
                 {'Ext Rot Pro',        'hip_rot',      'hip_flex',     0,      'knee_flex',    90}};
    LUT.K   =   {{'Flexion',            'knee_flex',    'hip_flex',     90}
                 {'Extension',          'knee_flex',    'hip_flex',     NaN}
                 {'Patella Alta'}
                 {'Popl Ang Uni',       'knee_flex',    'hip_flex',     90}
                 {'Popl Ang Bi',        'knee_flex',    'hip_flex',     90}
                 {'Spont Angle',        'knee_flex',    'hip_flex',     NaN}
                 {'Rectus Fem',         'knee_flex',    'hip_flex',     0}};
    LUT.A   =   {{'Dorsiflexion Kn 0',  'ankle_flex',   'knee_flex',    0}
                 {'Dorsiflexion Kn 90', 'ankle_flex',   'knee_flex',    90}
                 {'Plantarflexion',     'ankle_flex'}};
end
if isempty(Muscles)
    Muscles.H   =   {{'Flexion',            'glut_max1','glut_max2','glut_max3','glut_med1','glut_med2','glut_med3','glut_min1','glut_min2','glut_min3'}
                     {'Extension',          'iliacus','psoas'}
                     {'Abduction 0',        'add_mag2','add_mag3','add_long'}
                     {'Abduction 90',       'add_mag1','add_brev'}
                     {'Adduction',          ''}
                     {'Int Rot Sup',        'glut_med1','glut_med2','glut_med3','glut_min1','glut_min2','glut_min3'}
                     {'Int Rot Pro',        'glut_med1','glut_med2','glut_med3','glut_min1','glut_min2','glut_min3'}
                     {'Ext Rot Sup',        'glut_med1','glut_med2','glut_med3','glut_min1','glut_min2','glut_min3'}
                     {'Ext Rot Pro',        'glut_med1','glut_med2','glut_med3','glut_min1','glut_min2','glut_min3'}};
    Muscles.K   =   {{'Flexion',            'rectus_fem','vas_int','vas_med','vas_lat'}
                     {'Extension',          'semimem','semiten','sartorius','gracilis','bi_fem_lh'}
                     {'Patella Alta',       ''}
                     {'Popl Ang Uni',       'semimem','semiten','sartorius','gracilis','bi_fem_lh'}
                     {'Popl Ang Bi',        'semimem','semiten','sartorius','gracilis','bi_fem_lh'}
                     {'Spont Angle',        'semimem','semiten','sartorius','gracilis','bi_fem_lh'}
                     {'Rectus Fem',         'rectus_fem','vas_int','vas_med','vas_lat'}};
    Muscles.A   =   {{'Dorsiflexion Kn 0',  'soleus','gas_lat','gas_med'}
                     {'Dorsiflexion Kn 90', 'soleus','gas_lat','gas_med'}
                     {'Plantarflexion',     'tib_ant'}};
end
import org.opensim.modeling.*
MOD         =   Model(modelName);
SSS         =   MOD.initSystem();
SBE         =   MOD.getSimbodyEngine();
LRLR        =   {upper(Side)};%{'L','R'};
CS          =   MOD.getCoordinateSet();
CSs         =   CS.getSize();
MUSC        =   MOD.getMuscles();
Ms          =   MUSC.getSize();
MS          =   MOD.getMarkerSet();
AllMuscles  =   {};
for q=1:Ms
    AllMuscles{q,1}     =   MUSC.get(q-1).getName().toCharArray()';
end
% order defined according to the infromation from Pellenberg
%  0-    patella alta (this is actually to be performed with the GUI tool fro patella advancement...)
%  1-    hip extension
%  2-    hip flexion
%  3-    knee flexion (using hip extension)
%  4-    knee extension (using hip extension)
%  5-    knee spontaneous (using hip extension)
%  6-    rectus femoris (using hip extension)
%  7-    hip abduction 0 (using hip extension and knee extension)
%  8-    hip abduction 90 (using hip extension and knee flexion) < this has to be implemented in a more realistic way. For now just change the adduction value
%  9-    popliteal angle biarticular
% 10-    popliteal angle uniarticular
% 11-    ankle dorsiflexion 0 (using knee extension)
% 12-    ankle dorsiflexion 90
% 13-    hip internal rotation supine
% 14-    hip internal rotation prone (using hip extension)
% 15-    hip external rotation supine
% 16-    hip external rotation prone (using hip extension)
LMTout                  = [];
for lr=1%:2    
    LMT                 =   zeros(length(AllMuscles),1);
    %  1-    hip extension <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Extension');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipExtension(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); %#ok<*AGROW>
        % set the default hip extension value to be used in the following measurements
        defHext.(LRLR{lr})  =   max(-Minfo{2},0); % negative because the signs are opposite between clinic and opensim
    else
        % set the default hip extension value to be used in the following measurements
        defHext.(LRLR{lr})  =   0; % negative because the signs are opposite between clinic and opensim
    end
    %  2-    hip flexion <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Flexion');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipFlexion(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); %#ok<*AGROW>
        % set the default hip extension value to be used in the following measurements
        defHfle.(LRLR{lr})  =   min(-Minfo{2},90); % negative because the signs are opposite between clinic and opensim
    else
        % set the default hip extension value to be used in the following measurements
        defHfle.(LRLR{lr})  =   90; % negative because the signs are opposite between clinic and opensim
    end
    %  3-    knee flexion (using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Flexion');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneeFlexion(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); %#ok<*AGROW>
        % set the default hip extension value to be used in the following measurements
        defKfle.(LRLR{lr})  =   min(-Minfo{2},90); % negative because the signs are opposite between clinic and opensim
    else
        % set the default hip extension value to be used in the following measurements
        defKfle.(LRLR{lr})  =   90; % negative because the signs are opposite between clinic and opensim
    end
    %  4-    knee extension (using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Extension');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneeExtension(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); 
        % set the default hip extension value to be used in the following measurements
        defKext.(LRLR{lr})  =   max(-Minfo{2},0); % negative because the signs are opposite between clinic and opensim
    else
        % set the default hip extension value to be used in the following measurements
        defKext.(LRLR{lr})  =   0; % negative because the signs are opposite between clinic and opensim
    end
    %  5-    knee spontaneous (using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Spont Angle');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneeSpont(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    %  6-    rectus femoris (using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Rectus Fem');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneeRectFem(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    %  7-    hip abduction 0 (using hip extension and knee extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Abduction 0');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipAbduction0(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); %#ok<*AGROW>
    end
    %  8-    hip abduction 90 (using hip extension and knee flexion) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Abduction 90');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipAbduction90(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); %#ok<*AGROW>
    end
    %  9-    popliteal angle biarticular <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Popl Ang Bi');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneePAB(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 10-    popliteal angle uniarticular <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Knee.(LRLR{lr})(:,1),'Popl Ang Uni');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Knee.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.K,Measurement);
        Musc2use            =   strcat(Muscles.K{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   kneePAU(LUT.K{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 11-    ankle dorsiflexion 0 (using knee extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Ankle.(LRLR{lr})(:,1),'Dorsiflexion Kn 0');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Ankle.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.A,Measurement);
        Musc2use            =   strcat(Muscles.A{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   ankleDorsiflexion0(LUT.A{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 12-    ankle dorsiflexion 90 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Ankle.(LRLR{lr})(:,1),'Dorsiflexion Kn 90');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Ankle.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.A,Measurement);
        Musc2use            =   strcat(Muscles.A{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   ankleDorsiflexion90(LUT.A{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 13-    hip internal rotation supine <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Int Rot Sup');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipIntRotSup(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 14-    hip internal rotation prone(using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Int Rot Pro');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipIntRotPro(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 15-    hip external rotation supine <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Ext Rot Sup');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipExtRotSup(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % 16-    hip external rotation prone (using hip extension) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    I                   =   strcmpi(Values.Hip.(LRLR{lr})(:,1),'Ext Rot Pro');
    if any(I)
        setJointsZero();
        Minfo               =   Values.Hip.(LRLR{lr})(I,:);
        Measurement         =   Minfo{1};
        I2                  =   findIndexLUT(LUT.H,Measurement);
        Musc2use            =   strcat(Muscles.H{I2}(2:end),'_',lower(LRLR{lr}));
        Imusc               =   [];
        for q=1:length(Musc2use)
            Imusc           =   [Imusc find(strcmp(Musc2use{q},AllMuscles))];
        end
        LMT(Imusc,end+1)    =   hipExtRotPro(LUT.H{I2},lower(LRLR{lr}),Minfo,Musc2use); 
    end
    % pass the maximum values to the output
    LMTout              =   [LMTout max(LMT,[],2)];
end
LMTout                  =   max(LMTout,[],2);    

if USEPASSIVEMOMENTS
else
    PasMom=[];
end

%                               ~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~
%                               ~~~~~~~~~~~~|~~~~~~~~~~~~ other functions ~~~~~~~~~~~~|~~~~~~~~~~~~
%                               ~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~

    function I              =   findIndexLUT(lut,N)
        I       =   [];
        for qI=1:size(lut,1)
            if strcmpi(lut{qI}{1},N)
                I       =   qI;
                break
            end
        end
    end
    function []             =   setJointsZero()
        for c=1:CSs
            CS.get(c-1).setValue(SSS,0);
        end
    end
% functions for the hip
    function lmtHipFlex     =   hipFlexion(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi) 
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipFlex(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipExt      =   hipExtension(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi) % negative because the signs are opposite between clinic and opensim
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipExt(m,1)  =   M.getLength(SSS);
        end
    end
    function lmtHipAbd0     =   hipAbduction0(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs: kneeflex and hipflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defKext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defHext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipAbd0(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipAbd90    =   hipAbduction90(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs: kneeflex and hipflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defKfle.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipAbd90(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipAdd      =   hipAdduction(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipAdd(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipIRS      =   hipIntRotSup(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs: hipflex and kneeflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipIRS(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipIRP      =   hipIntRotPro(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs: hipflex and kneeflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defHext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipIRP(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipERS      =   hipExtRotSup(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs: hipflex and kneeflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipERS(m,1) =   M.getLength(SSS);
        end
    end
    function lmtHipERP      =   hipExtRotPro(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs: hipflex and kneeflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defHext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        for w=5
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtHipERP(m,1) =   M.getLength(SSS);
        end
    end
% functions for the knee
    function lmtKnFlex      =   kneeFlexion(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs (hip flexion)
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defHext.(upper(Side)))%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtKnFlex(m,1)  =   M.getLength(SSS);
        end
    end
    function lmtKnExt       =   kneeExtension(Info,Side,Minfo,Muscles)
        % the markers on the heel and psis are on the same line and the knee is flexed according to the value in the exam
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi) % value of knee extension is reported with an inverted sign
        % get information about the hip joint orientation and rotation axis
        Hdof        =   CS.get([Info{3},'_',Side]);
        HJ          =   Hdof.getJoint();
        OIP         =   org.opensim.modeling.Vec3(0);
        HJ.getOrientationInParent(OIP);
        OIP         =   [OIP.get(0) OIP.get(1) OIP.get(2)];
        OIP         =   Euler2mat(OIP);
        LIP         =   org.opensim.modeling.Vec3(0);
        HJ.getLocationInParent(LIP);
        Pelv        =   HJ.getParentBody();
        Ptra0       =   SBE.getTransform(SSS,Pelv).R;
        for r=1:3
            for c=1:3
                Ptra(r,c)   =   Ptra0.get(r-1,c-1);
            end
        end
        HJ          =   org.opensim.modeling.CustomJoint.safeDownCast(HJ);
        FlexST      =   HJ.getSpatialTransform().getTransformAxis(0).getAxis();
        FlexST      =   [FlexST.get(0) FlexST.get(1) FlexST.get(2)]';
        FlexAX      =   (OIP*FlexST);
        % get postition in body of the two markers to be used (on heel and pelvis)
        Mark1       =   MS.get([upper(Side),'HEE']);
        Mark2       =   MS.get([upper(Side),'PSIS']);
        M1posR      =   Mark1.getOffset();
        M2posR      =   Mark2.getOffset();
        M1Body      =   Mark1.getBody();
        M2Body      =   Mark2.getBody();
        M1pos0      =   org.opensim.modeling.Vec3(0);
        M2pos0      =   org.opensim.modeling.Vec3(0);
        LIP0        =   org.opensim.modeling.Vec3(0);
        % position of the markers relative to the pelvis reference frame
        SBE.getPosition(SSS,M1Body,M1posR,M1pos0)
        SBE.getPosition(SSS,M2Body,M2posR,M2pos0)
        SBE.getPosition(SSS,M2Body,LIP,LIP0)
        M1pos0      =   [M1pos0.get(0) M1pos0.get(1) M1pos0.get(2)]';
        M2pos0      =   [M2pos0.get(0) M2pos0.get(1) M2pos0.get(2)]';
        LIP0        =   [LIP0.get(0) LIP0.get(1) LIP0.get(2)]';
        % points in the pelvis reference frame
        LM2         =   Ptra'*(LIP0-M2pos0);
        M1LP        =   Ptra'*(M1pos0-LIP0);
        % define the angle to be applied to the hip, in order to bring the markers on the same lavel (X value)
        Hang        =   hipAngleClinicalExam(M1LP,FlexAX,-LM2(1));
        % check if the patient can reach the desired position
        Hang        =   max(Hang,defHext.(upper(Side))/180*pi);
        % change the postion of the hip to put heel and pelvis on the same plane
        Hdof.setValue(SSS,Hang)
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtKnExt(m,1)   =   M.getLength(SSS);
        end
    end
    function lmtKnSpo       =   kneeSpont(Info,Side,Minfo,Muscles)
        % the markers on the heel and psis are on the same line and the knee is flexed according to the value in the exam
    end
    function lmtKnPA        =   kneePatAlta(Info,Side,Minfo,Muscles)
        % to be defined
    end
    function lmtRecFemA     =   kneeRectFem(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        switch Minfo{2}
            case {'3',3}
                knAng   =   60;
            case {'2',2}
                knAng   =   90;
            case {'1',1}
                knAng   =   105;
            case {'0',0}
                knAng   =   120;
        end
        CS.get(nameDof).setValue(SSS,knAng/180*pi)
        % set the value of the other DoFs: hipflex
        for w=3
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defHext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtRecFemA(m,1) =   M.getLength(SSS);
        end
    end
    function lmtKnPAU       =   kneePAU(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtKnPAU(m,1)   =   M.getLength(SSS);
        end
    end
    function lmtKnPAB       =   kneePAB(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,-Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtKnPAB(m,1)   =   M.getLength(SSS);
        end
    end
% functions for the ankle 
    function lmtAnkDF0      =   ankleDorsiflexion0(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,defKext.(upper(Side))/180*pi)%CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtAnkDF0(m,1)  =   M.getLength(SSS);
        end
    end
    function lmtAnkDF90     =   ankleDorsiflexion90(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtAnkDF90(m,1) =   M.getLength(SSS);
        end
    end
    function lmtAnkPF       =   anklePlantarflexion(Info,Side,Minfo,Muscles)
        % seet the value of the main Dof
        nameDof     =   [Info{2},'_',Side];
        CS.get(nameDof).setValue(SSS,Minfo{2}/180*pi)
        % set the value of the other DoFs
        for w=3:2:length(Info)
            nameDof     =   [Info{w},'_',Side];
            CS.get(nameDof).setValue(SSS,Info{w+1}/180*pi)
        end
        % measure the lengths of the muscles
        for m=1:length(Muscles)
            M               =   MUSC.get(Muscles{m});
            lmtAnkPF(m,1)   =   M.getLength(SSS);
        end
    end

end