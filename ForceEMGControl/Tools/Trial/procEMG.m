function[trial_data] = procEMG(trial_data,Aparams)
% Initialize variables
fchEMG = 10;
fclEMG = [];
fnEMG = [];
fs = 1024;
avgWindow = 200;
EMGScale = [];

struct2vars(who,Aparams)

% Filtering parameters
wnh = (2/fs)*fchEMG;
[b,a] = butter(2,wnh,'high');

if ~isempty(fclEMG)
    wnl = (2/fs)*fclEMG;
    [d,c] = butter(2,wnl,'low');
end

if ~isempty(fnEMG)
    wnn = (2/fs)*fnEMG;
    [f,e] = iirnotch(wnn,wnn/35);
end

for i = 1:length(trial_data)
    EMGfiltH = filtfilt(b,a,trial_data(i).EMG.raw);
    
    if ~isempty(fnEMG)
        EMGfiltH = filtfilt(f,e,EMGfiltH);
    end

    EMGrect = abs(EMGfiltH);
    
    if ~isempty(fclEMG)
        EMGfiltL = filter(d,c,EMGrect);
    else
        EMGfiltL = EMGrect;
    end
    
    EMGavg = movingAvg(EMGfiltL,avgWindow);
    
    trial_data(i).EMG.filt = EMGfiltH;
    trial_data(i).EMG.rect = EMGfiltL;
    trial_data(i).EMG.avg = EMGavg;
    
    if isfield(trial_data(i).EMG,'offset') && (length(trial_data(i).EMG.offset) == size(trial_data(i).EMG.filt,2))
        trial_data(i).EMG.rectOffset = EMGfiltL./repmat(trial_data(i).EMG.offset,[size(trial_data(i).EMG.filt,1) 1]);
        trial_data(i).EMG.avgOffset = movingAvg(trial_data(i).EMG.rectOffset,avgWindow);
    end
end
end