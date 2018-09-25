function[MVCScale] = MVCtest(varargin)
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
    
    MVCScale = zeros(2,length(channelSubset)-1);
    
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
        samplesMVCRect = abs(samplesMVCFilt);
        samplesMVCFilt = filter(d,c,abs(samplesMVCFilt));
        
        MVCScale(1,j) = mean(samplesMVCFilt,2);
        MVCScale(2,j) = mean(samplesMVCRect,2);
    end
end

sampler.disconnect()
library.destroy()
end