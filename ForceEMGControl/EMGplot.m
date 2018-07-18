function[] = EMGplot(emg_data,sampler,channels,channel_subset)

emg_plot = TMSi.RealTimePlot('EMG RealTimePlot', sampler.sample_rate, channels);
emg_plot.setWindowSize(10);
emg_plot.show();

sampler.start();
samples = sampler.sample();
emg_data.append(samples(channel_subset,:));
while emg_plot.is_visible
    emg_plot.append(samples(channel_subset,:));
    emg_plot.draw();
end
sampler.stop();
end