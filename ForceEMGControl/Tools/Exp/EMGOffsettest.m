function[EMGOffset] = EMGOffsettest(varargin)
% EMG parameters
channelSubset =     [];
channelName =       {};
sampleRateEMG =     1024;
fchEMG =            20; % [Hz]
fclEMG =            60; % [Hz]

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

fprintf('\nRunning MVC calibration for EMG-control task.\n\n')

% Initializing EMG
library = TMSi.Library('usb');
[EMGEnabled,sampler,~,~] = EMGinit(library,channelSubset,channelName,sampleRateEMG);

if EMGEnabled
    fprintf('EMG initialized.\n\n')
    
    str = ['Press enter when prepared for EMG Offset calculation.'];
    input(str)
    
    EMGOffset = zeros(length(channelSubset)-1,2);
    samplesOffset = [];
    
    sampler.start()
    
    for n = 1:10
        samples = sampler.sample();
        pause(0.2)
        samplesOffset = [samplesOffset, samples(channelSubset(1:end-1),:)];
    end
    
    sampler.stop()
    
    wnh = (2/sampleRateEMG)*fchEMG;
    wnl = (2/sampleRateEMG)*fclEMG;
    [b,a] = butter(2,wnh,'high');
    [d,c] = butter(2,wnl,'low');
    samplesOffsetFilt = filtfilt(b,a,samplesOffset')';
    samplesOffsetFilt = filter(d,c,abs(samplesOffsetFilt),[],2);
    EMGOffset = mean(samplesOffsetFilt,2);
end
sampler.disconnect()
library.destroy()
end