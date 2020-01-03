% This function returns the IK given as inputs the path to the
% file containing the output of the IK and the joints of
% interest in the desired order

function IK = getIK_MRI(pathIK,joints)

IKall = importdata(pathIK);

IK.time = IKall.data(:,strcmp(IKall.colheaders,{'time'}));
IK.all(:,1) = IK.time;
IK.colheaders{1,1} = 'time';
count = 1;
for i = 1:length(joints)
    count = count + 1;
    if strcmp(joints{i},'lower_torso_TX') || strcmp(joints{i},'lower_torso_TY') || strcmp(joints{i},'lower_torso_TZ')
        IK.(joints{i}) = IKall.data(:,strcmp(IKall.colheaders,joints{i}));
    else
        IK.(joints{i}) = IKall.data(:,strcmp(IKall.colheaders,joints{i})).*(pi/180);
    end
    IK.all(:,count) = IK.(joints{i});
    IK.colheaders(1,count) = IKall.colheaders(1,strcmp(IKall.colheaders,joints{i}));
end
% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(IK.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
IK.allfilt = IK.all;
IK.allfilt(:,2:end) = filtfilt(af,bf,IK.allfilt(:,2:end));

end