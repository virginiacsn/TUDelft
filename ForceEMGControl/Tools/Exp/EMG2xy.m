function[EMGDatax,EMGDatay] = EMG2xy(emg_proc,ang)

if length(emg_proc) == 2
    xy = [cos(ang) -sin(ang); sin(ang) cos(ang)]*emg_proc;
    
    EMGDatax = xy(1);
    EMGDatay = xy(2);
else
    EMGDatax = [];
    EMGDatay = [];
end
end