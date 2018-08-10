function[EMGDatax,EMGDatay] = EMG2xy(emg_proc,ang)

if length(emg_proc) == 2
    xy = [cosd(ang) -sind(ang); sind(ang) cosd(ang)]*emg_proc;
    
    EMGDatax = xy(1);
    EMGDatay = xy(2);
else
    EMGDatax = [];
    EMGDatay = [];
end
end