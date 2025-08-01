function[trial_data] = fftEMG(trial_data)

for i = 1:length(trial_data)
    trial_data(i).EMG.fftraw = abs(fft(trial_data(i).EMG.raw))./size(trial_data(i).EMG.raw,1);
    trial_data(i).EMG.fftfilt = abs(fft(trial_data(i).EMG.filt))./size(trial_data(i).EMG.filt,1);
    trial_data(i).EMG.fftrect = abs(fft(trial_data(i).EMG.rect))./size(trial_data(i).EMG.rect,1);
end
end