
% This script will open the EMG and create and connect the sampler.
% parameters for sampling will be set
% Actual sampling is not done in this file. 
% That will be done in the Trial file
% Closing of the EMG device will also be done in the trial file

% Step 1: Setup library and choose connection type 
% Subset of channels to show
channel_subset = [1 2 3 4 17 18]; % only show channels 1,2,3,4

library = TMSi.Library('usb');

% Step 2: Find device to connect to. Keep on trying every second.
while numel(library.devices) == 0
    library.refreshDevices()
    pause(1);
end

try
    % Step 3: Get first device an retrieve information about device.
    device = library.getFirstDevice();
    
    % Step 4: Create a sampler with which we are going to retrieve samples.
    sampler = device.createSampler();
    
    
    
    try
        % Step 5: Set settings for sampler.
        sampler.setSampleRate(1024);
        %sampler.setBufferSize(0.5*SampleRate*device.channels)
        %data file?
        emg_data = TMSi.Data('Example', sampler.sample_rate, sampler.device.channels(channel_subset)); %subset of channels
        %emg_data = TMSi.Data('Example', sampler.sample_rate,sampler.device.channels);
        
        % sampler.setReferenceCalculation(true);
        
        % Step 7: Connect to device through the sampler.
        sampler.connect();
        EMGEnabled=1;
        messageEMG = 'EMG Initialized';
        
    catch matlabException
        warning(sprintf('Error while trying to sample.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
        EMGEnabled=0;
        messageEMG = 'EMG Initialization ERROR';
    end
    
catch matlabException
    warning(sprintf('Error while trying to create device/sampler.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
    EMGEnabled=0;
    messageEMG = 'EMG Initialization ERROR';
end

% %     % Step 11: Stop sampler.
% %     sampler.stop();
% % 
% %     % Step 12: Disconnect with device.
% %     sampler.disconnect();