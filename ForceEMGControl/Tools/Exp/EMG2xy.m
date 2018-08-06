function[EMGDatax,EMGDatay] = EMG2xy(emg_proc,numTargets)

if length(emg_proc) == 2
    if numTargets == 3
        xy = [cosd(45) -cosd(45); cosd(45) cosd(45)]*emg_proc;
    elseif numTargets == 2
        xy = [cosd(45) -cosd(45); cosd(45) -cosd(45)]*emg_proc;    
    end
elseif length(emg_proc) == 4
    xy = [cosd(45) -cosd(45) cosd(45) -cosd(45); cosd(45) cosd(45) -cosd(45) -cosd(45)]*emg_proc;
else
    xy = [0 0];
end

EMGDatax = xy(1);
EMGDatay = xy(2);

end