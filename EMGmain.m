
channel_subset = [1 3 18];
channel_name = {'BB','TB','Saw'};
if length(channel_subset)~=length(channel_name)
    error('Names for all channels not available.')
end
sample_rate = 1024;
plotEMG = 0;

EMGinit;

if EMGEnabled
    sampler.start()
    
    if plotEMG
        emg_plot = TMSi.RealTimePlot('RealTimePlot Example', sampler.sample_rate, channels);
        emg_plot.setWindowSize(10);
        emg_plot.show();
    end
    
    tstart = tic;
    
    while toc(tstart)<5
        samples = sampler.sample();
        emg_data.append(samples(channel_subset,:));
        if emg_plot.is_visible
            emg_plot.append(samples(channel_subset,:));
            emg_plot.draw();
        end
    end
    
    EMGdataout = emg_data.samples;
    sampler.stop()
    sampler.disconnect()
    
    library.destroy()
    
else
    EMGdataout = [];
end
