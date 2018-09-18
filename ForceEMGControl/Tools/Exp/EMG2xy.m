function[EMGDatax,EMGDatay] = EMG2xy(emg_proc,ang)

if length(ang) == 1
    xy = [cos(ang) -sin(ang); sin(ang) cos(ang)]*emg_proc;
    
    EMGDatax = xy(1);
    EMGDatay = xy(2);
elseif length(ang) == 4
    
else
    EMGDatax = 0;
    EMGDatay = 0;
end
end