%% DAQ from EMG - CO Task
% Virginia Casasnovas
% 11/07/2018

function EMGControl_CO(varargin)
%% Parameter assignment
% Default parameters
% File parameters
saveforce =  0;
saveEMG =    1;
date =      '20180717';
task =      'CO';
code =      '001';
filenameforce = [date,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filepath =  [pwd '\Data\' date '\'];

% Task parameters
numTrials =         3;

targetForce =       1000; % [N]
numTargets =        3;
rCirTarget =        targetForce/10; % [N]
rCirCursor =        targetForce/20; % [N]

fchEMG =            10; % [Hz]

scanRate =          1000; % [scans/sec]
availSamples =      100; % [samples]
bufferWin =         200; % [samples]
iterUpdatePlot =    10;

movemtime =         5; % sec
holdtime =          1; % sec
timeout =           1; % sec
relaxtime =         1; % sec

% EMG parameters
plotEMG =           0;
channelSubset =     [1 2];
channelName =       {'BB','TL'};
sampleRateEMG =     1024;
smoothWin =         500;

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

if length(channelSubset)~=length(channelName)
    error('Names for all channels not available.')
end

%% Initialization
disp('Running DataAcquisition for CO EMG task.')

library = TMSi.Library('usb');
[EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channelSubset,channelName,sampleRateEMG);

if EMGEnabled
    disp('EMG initialized.')
    
    if saveEMG || saveforce
        if exist([filepath,filenameEMG],'file')||exist([filepath,filenameforce],'file')
            savefile = input('\nThis file already exsists. Continue saving? (y/n) ','s');
            if strcmp(savefile,'y')
                if saveEMG
                    fprintf('Saving EMG data in %s.\n',filenameEMG)
                end
                if saveforce
                    fprintf('Saving force data in %s.\n\n',filenameforce)
                else
                    fprintf('\n')
                end
            else
                saveEMG = 0;
                saveforce = 0;
                device = [];
                fprintf('Not saving data.\n\n')
            end
        else
            if saveEMG
                fprintf('Saving EMG data in %s.\n',filenameEMG)
            end
            if saveforce
                fprintf('Saving force data in %s.\n\n',filenameforce)
            else
                fprintf('\n')
            end
        end
    else
        fprintf('Not saving data.\n\n')
    end
    
    % Get EMG offset
    input('Press enter when prepared for EMG offset calculation')
    samplesOffset = [];
    sampler.start()
    for n = 1:20
        samples = sampler.sample();
        samplesOffset = [samplesOffset, samples(channelSubset,:)];
    end
    sampler.stop()
    EMGOffset = mean(samplesOffset,2);
    
    MVCScale = zeros(length(channelSubset),1);
    for j = 1:length(channelSubset)
        str = ['Press enter when prepared for ',channelName{j},' EMG MVC calculation'];
        input(str)
        samplesMVC = [];
        sampler.start()
        for n = 1:20
            samples = sampler.sample();
            samplesMVC = [samplesMVC, samples(channelSubset(j),:)];
        end
        sampler.stop()
        MVCScale(j) = mean(samplesMVC);
    end
    fprintf('\n')
    
    % Plot EMG
    if plotEMG
        disp('Plotting EMGs to check.')
        EMGplot(emg_data,sampler,channels,channelSubset); %emg_data appended will be saved? trigger
    end
    
    % Initialize force DAQ
    if saveforce
        
        device = daq.getDevices;
        
        if ~isempty(device)
            % Create NI DAQ session
            disp('Creating NI DAQ session.')
            s = daq.createSession('ni');
            addAnalogInputChannel(s,device.ID,0:6,'Voltage');
            addAnalogOutputChannel(s,device.ID,'ao0','Voltage');
            
            % Obtain offset by averaging 2 sec of still data
            fprintf('\nObtaining offset values...\n')
            queueOutputData(s,zeros(2*scanRate,1));
            forceOffsetData = s.startForeground();
            forceOffset = mean(forceOffsetData(:,1:6),1);
            disp('Offset obtained.')
            
            % Save force file header 
            if saveforce
                samplenum = 1;
                forceDataOut(samplenum,:) = {'Trialnum', 'State', 'TimeStamp', 'Fx', 'Fy', 'Fz','Trigger'};
            end
        else
            device = [];
            disp('DAQ device not found.')
        end
    end
    
    %% Data acquisition
    input('\nPress enter to start acquisition.')
    
    % Initialize variables
    global htrl
    
    countState = 0;
    
    % Set target forces
    targetAngles = [pi/4:pi/(2*(numTargets-1)):3*pi/4]; % [rad]
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    
    % Set figure
    hf = figure('Name','CO EMG Control Task');
    [hf,hp] = Figinit(hf,targetForce);
    title('EMG');
    xlabel(channelName{1}); ylabel(channelName{2});
    xl = xlim; yl = ylim;
    htrl = text(xl(2)+0.3*xl(2),yl(2),'Trial: 0','clipping','off','Fontsize',14);
    %set(gca,'XTickLabel',[],'YTickLabel',[])
    
    % Start EMG data sampling
    sampler.start()
    
    if ~isempty(device)
        % Initialize variables
        global tmove trelax tfail tsuccess tholdstart
        global targetCir
        global htrg hsta
        
        trialNum = 0;
        EMGDataBuffer = zeros(length(channelSubset),smoothWin);
        countBuffer = 0;
        cursorHoldOut = 0;
        state = 'start';
        tempState = 'start';
        emg_save = [];
        
        % Save EMG file header
        if saveEMG
            samplenum = 1;
            EMGDataOut(samplenum,:) = {'Trialnum', 'State', 'EMG'};
        end
        
        % Add event listener and start acquisition
        outputData = [4*ones(1,availSamples),zeros(1,scanRate-availSamples)]';
        queueOutputData(s,outputData);
        hlout = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));
        
        hlin = addlistener(s,'DataAvailable',@(src,event) processForceData(event,forceOffset,EMGOffset,MVCScale,hp));
        s.IsContinuous = true;
        s.Rate = scanRate; % scans/sec, samples/sec?
        s.NotifyWhenDataAvailableExceeds = availSamples; % Call listener when x samples are available
        s.startBackground();

        input('\Press enter to stop acquisition.') 

    else
        disp('No DAQ device detected.')
        
        for itrial = 1:numTrials
            
            nextTrial = false;
            cursorHoldOut = 0;
            state = 'start';
            tempState = 'start';
            
            while ~nextTrial
                
                samples = sampler.sample();
                emg_data.append(samples(channelSubset,:));

                tvis = tic;
                pause(0.02)
                
                if toc(tvis) >= 0.02
                    emg_data_len = size(emg_data.samples,2);
                    
                    % fprintf('%d\n',size(samples_subset,2));
        
                    wn = (2/sampleRateEMG)*fchEMG;
                    [b,a] = butter(2,wn,'high');
                    filtemg = filter(b,a,emg_data.samples-EMGOffset);
                    emg_proc = mean(abs(filtemg),2); % Rectify and average
                    
                    [EMGDatax,EMGDatay] = EMG2xy(emg_proc);
                    
                    cursorCir = circle(rCirCursor,EMGDatax,EMGDatay);
                    set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
                    
                    if strcmp(state,tempState) && countState == 0
                        countState = countState+1;
                        xl = xlim;
                        hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
                    elseif ~strcmp(state,tempState) && countState > 0
                        countState = 0;
                        delete(hsta)
                    end
                    
                    drawnow;
                    
                    % Trial
                    tempState = state;
                    
                    switch state
                        case 'start'
                            xl = xlim; yl = ylim;
                            delete(htrl)
                            htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(itrial)],'clipping','off','Fontsize',14);
                            
                            iAngle = randi(numTargets);
                            targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                            htrg = plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
                            
                            state = 'movement';
                            tmove = tic;
                        case 'movement'
                            if cursorInTarget(cursorCir,targetCir) && toc(tmove) <= movemtime
                                state = 'hold';
                                tholdstart = tic;
                            elseif ~cursorInTarget(cursorCir,targetCir) && toc(tmove) > movemtime
                                state = 'fail';
                                tfail = tic;
                            end
                        case 'hold'
                            if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart) <= holdtime
                                cursorHoldOut = cursorHoldOut+1;
                            elseif cursorHoldOut > 5 && toc(tholdstart) > holdtime
                                cursorHoldOut = 0;
                                state = 'fail';
                                tfail = tic;
                            elseif cursorHoldOut <= 5 && toc(tholdstart) > holdtime
                                cursorHoldOut = 0;
                                state = 'success';
                                tsuccess = tic;
                            end
                        case 'fail'
                            if toc(tfail) > timeout
                                state = 'relax';
                                trelax = tic;
                                delete(htrg)
                            end
                        case 'success'
                            if toc(tsuccess) > timeout
                                state = 'relax';
                                trelax = tic;
                                delete(htrg)
                            end
                        case 'relax'
                            if cursorInTarget(cursorCir,circle(1.5*rCirCursor,0,0)) && toc(trelax) > relaxtime
                                state = 'start';
                                nextTrial = true;
                            end
                    end
                end
            end
        end
        
        % Appending trial data
        if saveEMG
            samplenum = samplenum+1;
            EMGDataOut(samplenum,:) = {itrial,state,samples_subset'};
        end
    end
    
    % Delete handles and stop session
    close(hf)
    if ~isempty(device)
        s.stop()
        delete(hlin)
        delete(hlout)
    end
    
    % Close EMG
    sampler.stop()
    sampler.disconnect()
    
    % Save data
    if saveforce && ~isempty(device)
        save([filepath,filenameforce],'forceDataOut')
    end
    if saveEMG
        EMGDataOut = emg_data.samples;
        save([filepath,filenameEMG],'EMGDataOut')
        %save('emg_proc','emg_save')
    end
    
else
    error('EMG could not be initialized.')
end

library.destroy()

%% Function handles
    function processForceData(event,forceOffset,EMGOffset,MVCScale,hp)
        calibMat = [-0.56571    -0.01516    1.69417     -31.81016   -0.35339    33.73195;...
            -0.67022    37.27575    1.14848     -18.75980   0.34816     -19.33789;...
            19.06483    0.45530     18.57588    0.72627     19.51260    0.65624;...
            0.71302     0.12701     -32.12926   -1.43043    32.51061    1.09333;...
            37.80563    1.05190     -19.89516   -0.79823    -18.79634   -0.79836;...
            0.48287     -18.59045   0.38579     -18.05502   0.75146     -19.86235];
        offsetMat = repmat(forceOffset,size(event.Data(:,1:6),1),1);
        offsetData = event.Data(:,1:6)-offsetMat;
        procData = zeros(size(event.Data(:,1:6)));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        timeStamp = event.TimeStamps;
        triggerData = event.Data(:,7);
        forceData = procData(:,1:3);
        
        samples = sampler.sample();
        nSamples = size(samples,2);
        appendSamples = (samples(channelSubset,:)-repmat(EMGOffset,1,nSamples))./repmat(MVCScale,1,nSamples);
        
        emg_data.append(appendSamples)
        
        bufferTemp = EMGDataBuffer;
        bufferTemp(:,1:nSamples) = (samples(channelSubset,:)-repmat(EMGOffset,1,nSamples))./repmat(MVCScale,1,nSamples);
        bufferTemp(:,nSamples+1:end) = EMGDataBuffer(:,1:smoothWin-nSamples);
        EMGDataBuffer = bufferTemp;
        countBuffer = countBuffer+1;
       
        wn = (2/sampleRateEMG)*fchEMG;
        [b,a] = butter(2,wn,'high');
        filtEMG = filter(b,a,EMGDataBuffer,[],2);
        emg_proc = mean(abs(filtEMG),2); % Rectify and average
        %emg_save = [emg_save,emg_proc];
        [EMGDatax,EMGDatay] = EMG2xy(emg_proc);
        
        cursorCir = circle(rCirCursor,EMGDatax,EMGDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if strcmp(state,tempState) && countState == 0
            countState = countState+1;
            xl = xlim;
            hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
        elseif ~strcmp(state,tempState) && countState > 0
            countState = 0;
            delete(hsta)
        end
        
        if countBuffer == iterUpdatePlot
            drawnow;
            countBuffer = 0;
        end
        
        tempState = state;
        
        % Trial
        switch state
            case 'start'                
                xl = xlim; yl = ylim;
                delete(htrl)
                htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(trialNum)],'clipping','off','Fontsize',14);
                
                iAngle = randi(numTargets);
                targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                htrg = plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
                
                state = 'movement';
                tmove = tic;
            case 'movement'
                if cursorInTarget(cursorCir,targetCir) && toc(tmove) <= movemtime
                    state = 'hold';
                    tholdstart = tic;
                elseif ~cursorInTarget(cursorCir,targetCir) && toc(tmove) > movemtime
                    state = 'fail';
                    tfail = tic;
                end
            case 'hold'
                if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart) <= holdtime
                    cursorHoldOut = cursorHoldOut+1;
                elseif cursorHoldOut > 5 && toc(tholdstart) > holdtime
                    cursorHoldOut = 0;
                    state = 'fail';
                    tfail = tic;
                elseif cursorHoldOut <= 5 && toc(tholdstart) > holdtime
                    cursorHoldOut = 0;
                    state = 'success';
                    tsuccess = tic;
                end
            case 'fail'
                if toc(tfail) > timeout
                    state = 'relax';
                    trelax = tic;
                    delete(htrg)
                end
            case 'success'
                if toc(tsuccess) > timeout
                    state = 'relax';
                    trelax = tic;
                    delete(htrg)
                end
            case 'relax'
                if cursorInTarget(cursorCir,circle(1.5*rCirCursor,0,0)) && toc(trelax) > relaxtime
                    state = 'start';
                    trialNum = trialNum + 1;
                end
        end
                
        % Appending trial data
        if saveforce
            samplenum = samplenum+1;
            forceDataOut(samplenum,:) = {trialNum,state,timeStamp,forceData(:,1),forceData(:,2),forceData(:,3),triggerData};
        end
    end
end