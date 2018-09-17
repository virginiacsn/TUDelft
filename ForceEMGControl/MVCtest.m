function[MVCScale] = MVCtest(varargin)
% EMG parameters
channelSubset =     [1 2 17];
channelControl =    [1 2];
channelName =       {'BB','TL','Trigger'};
sampleRateEMG =     1024;
fchEMG =            20; % [Hz]
fclEMG =            60;

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

fprintf('\nRunning MVC calibration for EMG-control task.\n\n')

% Initializing EMG
library = TMSi.Library('usb');
[EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channelSubset,channelName,sampleRateEMG);

if EMGEnabled
    fprintf('EMG initialized.\n\n')
    
    MVCScale = zeros(length(channelSubset)-1,1);
    
    for j = 1:length(channelSubset)-1
        str = ['Press enter when prepared for ',channelName{channelSubset(j)},' EMG MVC calculation.'];
        input(str)
        samplesMVC = [];
        
        sampler.start()
        
        for n = 1:10
            samples = sampler.sample();
            pause(0.2)
            samplesMVC = [samplesMVC, samples(channelSubset(j),:)];
        end
        
        sampler.stop()
        
        wnh = (2/sampleRateEMG)*fchEMG;
        wnl = (2/sampleRateEMG)*fclEMG;
        [b,a] = butter(2,wnh,'high');
        [d,c] = butter(2,wnl,'low');
        
        samplesMVCFilt = filtfilt(b,a,samplesMVC);
        samplesMVCFilt = filter(d,c,samplesMVCFilt);
        MVCScale(j) = mean(abs(samplesMVCFilt),2);
    end
end
sampler.disconnect()
library.destroy()
end