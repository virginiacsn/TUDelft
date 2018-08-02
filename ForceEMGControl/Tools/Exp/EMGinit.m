function[varargout] = EMGinit(library,channel_subset,channel_name,sample_rate)
% Step 2: Find device to connect to. Keep on trying every second.
tscan = tic;
while numel(library.devices) == 0 && toc(tscan) < 10
    library.refreshDevices();
    pause(1);
end

if numel(library.devices) == 0
    disp('EMG device not found.')
    EMGEnabled = 0;
    varargout = {EMGEnabled,[],[],[]};
else
    
    try
        % Step 3: Get first device an retrieve information about device.
        device = library.getFirstDevice();
        
        % Step 4: Create a sampler with which we are going to retrieve samples.
        sampler = device.createSampler();
        
        try
            % Step 5: Set settings for sampler.
            sampler.setSampleRate(sample_rate);
            % sampler.setReferenceCalculation(true);
            channels = device.channels(channel_subset);
            
            if ~isempty(channel_name)
                for i = 1:length(channel_subset)
                    channels{i}.name = channel_name{i};
                end
            end
            
            emg_data = TMSi.Data('Example', sampler.sample_rate, sampler.device.channels(channel_subset)); %subset of channels
            
            % Step 7: Connect to device through the sampler.
            sampler.connect();
            EMGEnabled = 1;
            varargout = {EMGEnabled,sampler,emg_data,channels};
            
            %         % Step 8: Start sampling.
            %         sampler.start();
            %
            %         % Step 9: Sample as long as plot is visible.
            % 				%   Can be closed with 'q' or window close cross.
            % 		        %   All ranges can be set by 'r', a dialog allows to put in a range.
            % 				%   Press 'a' for autoscale.
            %         while pltEMG.is_visible
            %             samples = sampler.sample();
            %             % Step 10: Append samples to plot.
            %             pltEMG.append(samples);
            %             pltEMG.draw();
            %         end
            
        catch matlabException
            warning(sprintf('Error while trying to sample.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
        end
        
        %     % Step 11: Stop sampler.
        %     sampler.stop();
        %
        %     % Step 12: Disconnect with device.
        %     sampler.disconnect();
        
    catch matlabException
        warning(sprintf('Error while trying to create device/sampler.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
        EMGEnabled = 0;
        varargout = {EMGEnabled,[],[],[]};
    end
end
end

% % Step 13: Cleanup library.
% library.destroy();
