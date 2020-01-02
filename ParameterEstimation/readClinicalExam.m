function [varargout]=readClinicalExam(varargin)
fname                   =   varargin{1};
Side                    =   varargin{2};
% define the correspondences between the names in the GUI and in the excel filie
LUT.H                   =   {'Flexion',             'pROM_Hpfl_',            'ROM_Hpfl'
                             'Extension',           'pROM_Hpext_supine_',    'ROM_Hpext_supine'
                             'Abduction 0',         'p_ROM_Hpabd0_',         'ROM_Hpabd0'
                             'Abduction 90',        'p_ROM_Hpabd90_',        'ROM_Hpabd90'
                             'Adduction',           'p_ROM_Hpadd_',          'ROM_Hpadd'
                             'Int Rot Sup',         'p_ROM_Hpendosup_',      'ROM_Hpendosup'
                             'Int Rot Pro',         'p_ROM_Hpendopro_',      'ROM_Hpendopro'
                             'Ext Rot Sup',         'p_ROM_Hpexosup_',       'ROM_Hpexosup'
                             'Ext Rot Pro',         'p_ROM_Hpexopro_',       'ROM_Hpexopro'      };
LUT.K                   =   {'Flexion',             'p_ROM_Knfl_',           'ROM_Knfl'
                             'Extension',           'p_ROM_Knext_',          'ROM_Knext'
                             'Popl Ang Uni',        'p_ROM_Popluni_',        'ROM_Popluni'
                             'Popl Ang Bi',         'p_ROM_Poplbi_',         'ROM_Poplbi'
                             'Spont Angle',         'p_ROM_Knsponpos_',      'ROM_Knsponpos'
                             'Rectus Fem',          'p_ROM_Ref_',            'ROM_Ref'           };
LUT.A                   =   {'Dorsiflexion Kn 0',   'pROM_Ankledf 0_',       'ROM_Ankledf 0'
                             'Dorsiflexion Kn 90',  'pROM_Ankledf90_',       'ROM_Ankledf90'     };
[~,~,tmpD]              =   xlsread(fname);
Head                    =   tmpD(1,:);
Data                    =   tmpD(2:end,:);
% find the date of the exam
indDate                 =   strcmp(Head,'MD');
dates                   =   Data(:,indDate);
if length(dates)>1
    [dateSel,~]         =   listdlg('PromptString','Select the date to be used:','SelectionMode','single','ListString',dates);
    dates               =   dates(dateSel); %#ok<NASGU>
    Data                =   Data(dateSel,:);
end
% empty structure for the values
readVal.Hip.L           =   {};
readVal.Hip.R           =   {};
readVal.Knee.L          =   {};
readVal.Knee.R          =   {};
readVal.Ankle.L         =   {};
readVal.Ankle.R         =   {};
% get the values from the excel and put them into the format for the GUI
for q=1:size(LUT.H,1)
    ind1                =   ~cellfun(@isempty,regexpi(Head,[LUT.H{q,3},'_',Side]));
    ind2                =   strncmpi(Head,'p',1);
    ind                 =   ind1&ind2;
    if sum(ind)>1
        warning(['More than one value for ',[LUT.H{q,1},' ',Side],' -> ',[LUT.H{q,2},Side]])
        tmp             =   find(ind,1,'first');
        ind             =   false(size(ind));
        ind(tmp)        =   true;
    end
    if any(ind)
        if ~isempty(Data{ind})&&isa(Data{ind},'double')&&~isnan(Data{ind})
            readVal.Hip.(Side){end+1,1} =   LUT.H{q,1};
            readVal.Hip.(Side){end,2}   =   Data{ind};
            readVal.Hip.(Side){end,3}   =   0;
        end
    end
end
for q=1:size(LUT.K,1)
    ind1                =   ~cellfun(@isempty,regexpi(Head,[LUT.K{q,3},'_',Side]));
    ind2                =   strncmpi(Head,'p',1);
    ind                 =   ind1&ind2;
    if sum(ind)>1
        warning(['More than one value for ',[LUT.K{q,1},' ',Side],' -> ',[LUT.K{q,2},Side]])
        tmp             =   find(ind,1,'first');
        ind             =   false(size(ind));
        ind(tmp)        =   true;
    end
    if any(ind)
        if ~isempty(Data{ind})&&isa(Data{ind},'double')&&~isnan(Data{ind})
            readVal.Knee.(Side){end+1,1} =   LUT.K{q,1};
            readVal.Knee.(Side){end,2}   =   Data{ind};
            readVal.Knee.(Side){end,3}   =   0;
        end
    end
end
for q=1:size(LUT.A,1)
    ind1                =   ~cellfun(@isempty,regexpi(Head,[LUT.A{q,3},'_',Side]));
    ind2                =   strncmpi(Head,'p',1);
    ind                 =   ind1&ind2;
    if sum(ind)>1
        warning(['More than one value for ',[LUT.A{q,1},' ',Side],' -> ',[LUT.A{q,2},Side]])
        tmp             =   find(ind,1,'first');
        ind             =   false(size(ind));
        ind(tmp)        =   true;
    end
    if any(ind)
        if ~isempty(Data{ind})&&isa(Data{ind},'double')&&~isnan(Data{ind})
            readVal.Ankle.(Side){end+1,1} =   LUT.A{q,1};
            readVal.Ankle.(Side){end,2}   =   Data{ind};
            readVal.Ankle.(Side){end,3}   =   0;
        end
    end
end

varargout{1}            =   readVal;


