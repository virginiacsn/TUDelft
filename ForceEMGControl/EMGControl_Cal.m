% Virginia Casasnovas
% 18/09/2018

% Function for calibration using EMG-control task. Center-out type. Saves
% EMG and force data in separate txt files.
% Input: structs with parameters: fileparams, taskparams, EMGparams,
% forceparams. Input values will overwrite defaults.

function EMGControl_Cal(varargin)
%% Parameter assignment
% Default parameters
% File parameters
saveforce =         0;
saveEMG =           0;
filenameforce =     [];
filenameEMG =       [];
filepath =          [];

% Task parameters
targetEMGCal =      50; 
targetTolEMG =      0.2;
cursorTolEMG =      4;

movemtimeCal =      2; % [sec]
holdtimeCal =       4; % [sec]
timeout =           1; % [sec]
relaxtime =         1; % [sec]

setFig =            1;

% Force parameters
scanRate =          1000; % [scans/sec]
availSamplesEMG =   500; % [samples]

% EMG parameters
plotEMG =           0;
channelSubset =     [];
channelSubsetCal =  [];
channelName =       {};
channelNameCal =    {};
channelAngleCal =   [];
sampleRateEMG =     1024;
fchEMG =            10; % [Hz]
fclEMG =            60;
smoothWin =         500;
iterUpdatePlotEMG = 1;
EMGScale =          [];

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

rCirTarget = targetEMGCal*targetTolEMG; % [N]
rCirCursor = targetEMGCal*targetTolEMG/cursorTolEMG; % [N]

if length(channelSubset)~=length(channelName)
    error('Names for all channels not available.')
end

if length(channelSubsetCal)~=length(channelNameCal)
    error('Names for all channels not available.')
end

[targetAnglesEMG,isort] = sort(channelAngleCal);
channelSubsetTemp = channelSubsetCal(1:end-1);
channelSubsetCal = [channelSubsetTemp(isort) channelSubsetCal(end)];
channelSubsetTemp = [channelSubsetTemp channelSubsetCal(end)];
% channelNameTemp = channelNameCal(1:end-1);
% channelNameCal = [channelNameTemp(isort) channelNameCal(end)];

%% Saving file check
if saveEMG || saveforce
    if exist([filepath,filenameEMG],'file')||exist([filepath,filenameforce],'file')
        savefile = input(['\n',filenameEMG,' already exsists. Continue saving? (y/n) '],'s');
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

%% Initialization
disp('Running DataAcquisition for calibration with EMG task.')

library = TMSi.Library('usb');
[EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channelSubset,channelName,sampleRateEMG);

if EMGEnabled
    disp('EMG initialized.')
    
    % Plot EMG
    if plotEMG
        disp('Plotting EMGs to check.')
        EMGplot(emg_data,sampler,channels,channelSubsetCal); %emg_data appended will be saved? trigger
    end
    
    % Initialize force DAQ
    if saveforce
        device = daq.getDevices;
        if length(device) > 1
            device = device(2);
        end
        if ~isempty(device)
            % Create NI DAQ session
            disp('Creating NI DAQ session.')
            s = daq.createSession('ni');
            addAnalogInputChannel(s,device.ID,0:6,'Voltage');
            addAnalogOutputChannel(s,device.ID,'ao0','Voltage');
            
            % Obtain offset by averaging 2 sec of still data
            input('Press enter when prepared for sensor offset calculation.')
            fprintf('Obtaining offset values...\n')
            queueOutputData(s,zeros(2*scanRate,1));
            forceOffsetData = s.startForeground();
            forceOffset = mean(forceOffsetData(:,1:6),1);
            disp('Offset obtained.')
            
            % Save force file header
            sampleNum = 1;
            forceDataOut_EMGCO(sampleNum,:) = {'Trialnum', 'TargetAng', 'State', 'TimeStamp', 'Fx', 'Fy', 'Fz','Trigger'};
            EMGDataOut_EMGCO(sampleNum,:) = {'EMG'};
        else
            device = [];
            disp('DAQ device not found.')
        end
    else
        device = [];
        disp('Not saving force data.')
    end
    
    % Obtain EMG offset
    repeat = 'y';
    while strcmp(repeat,'y')
        input('\nPress enter when prepared for EMG offset calculation.')
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
        
        EMGOffsetCal = EMGOffset(channelSubsetCal(1:end-1));
        
        fprintf('EMG offset:\n')
        for k = 1:length(channelSubsetCal)-1
            fprintf('%s: %1.3f\n',channelNameCal{channelSubsetTemp==channelSubsetCal(k)},EMGOffsetCal(k))
        end
        fprintf('\n')
        repeat = input('Repeat EMG offset calculation? [y/n] ','s');
    end
    fprintf('\n')
    
    %% Data acquisition
    input('\nPress enter to start acquisition.')
    
    % Initialize variables
    global htrl 
    
    countState = 0;
    
    % Set target forces
    targetPosx = targetEMGCal*cos(targetAnglesEMG)';
    targetPosy = targetEMGCal*sin(targetAnglesEMG)';
    
    % Set figure
    hf = figure('Name','CO EMG Control Task');
    %[hf,hp] = Figinit(hf,[max(targetPosx) max(targetPosy)]./1.2);
    [hf,hp] = Figinit(hf,targetEMGCal*[1 1]);
    if setFig
        title('EMG');
        %xlabel(channelName{channelControl(1)}); ylabel(channelName{channelControl(2)});
    else
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        set(gca,'Color','k');
    end
    xl = xlim; yl = ylim;
    htrl = text(xl(2)+0.3*xl(2),yl(2),'Trial: 0','clipping','off','Fontsize',14);
    
    % Start EMG data sampling
    sampler.start()
    
    if ~isempty(device)
        % Initialize variables
        global tmove trelax tfail tsuccess tholdstart
        global targetCir iAngle trialNum
        global htrg hsta 
        
        trialNum = 0;
        countBuffer = 0;
        EMGDataBuffer = zeros(length(channelSubsetCal)-1,smoothWin);
        state = 'start';
        tempState = 'start';
        cursorHoldOut = 0;
        iAngle = 1;
        reptarg = 0;
        
        emg_save = [];
        
        % Add event listener and start acquisition
        hlin = addlistener(s,'DataAvailable',@(src,event) processForceData(event,forceOffset,EMGOffsetCal,EMGScale,hp));
        %hstop = addlistener(s,'DataAvailable',@(src,event) stopTrialNum(trialNum,numTrials));
        s.IsContinuous = true;
        s.Rate = scanRate; % scans/sec, samples/sec?
        s.NotifyWhenDataAvailableExceeds = availSamplesEMG; % Call listener when x samples are available
        
        outputData = [4*ones(1,availSamplesEMG),zeros(1,scanRate-availSamplesEMG)]';
        queueOutputData(s,outputData);
        hlout = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));
        
        s.startBackground();
        
        input('\Press enter to stop acquisition.');
    
    end
    
    % Delete handles and stop session
    close(hf)
    if ~isempty(device)
        s.stop()
        delete(hlin)
        delete(hlout)
        %delete(hstop)
    end
    
    % Close EMG
    sampler.stop()
    sampler.disconnect()
    
    % Save data
    if ~isempty(device)
        if saveforce
            save([filepath,filenameforce],'forceDataOut_EMGCO')
        end
    end
    if saveEMG
        EMGOffset = EMGOffset';
        save([filepath,filenameEMG],'EMGDataOut_EMGCO','EMGOffset','EMGScale')
    end
    
else
    error('EMG could not be initialized.')
end

library.destroy()

%% Function handles
    function processForceData(event,forceOffset,EMGOffsetCal,EMGScale,hp)
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
        appendSamples = samples(channelSubset,:);%-repmat(EMGOffset,1,nSamples))./repmat(EMGScale,1,nSamples);
        
        %emg_data.append(appendSamples)
        EMGSamples = samples(channelSubsetCal(1:end-1),:);
        
        if nSamples < smoothWin
            bufferTemp = EMGDataBuffer;
            bufferTemp(:,1:nSamples) = EMGSamples;
            bufferTemp(:,nSamples+1:end) = EMGDataBuffer(:,1:smoothWin-nSamples);
            EMGDataBuffer = bufferTemp;
            countBuffer = countBuffer+1;
        end
        
        wnh = (2/sampleRateEMG)*fchEMG;
        wnl = (2/sampleRateEMG)*fclEMG;
        [b,a] = butter(2,wnh,'high');
        [d,c] = butter(2,wnl,'low');
        
        filtEMGBuffer = filtfilt(b,a,EMGDataBuffer')';
        filtEMGBuffer = filter(d,c,abs(filtEMGBuffer),[],2);
 
        avgRectEMGBuffer = (mean((filtEMGBuffer),2)-EMGOffsetCal)./(EMGScale(channelSubsetCal(1:end-1))-EMGOffsetCal); % Rectify, smooth and scale
        avgRectEMGBuffer(isnan(avgRectEMGBuffer)) = 0;
        emg_save = [emg_save,avgRectEMGBuffer];
        
        [EMGDatax,EMGDatay] = EMG2xy([avgRectEMGBuffer(iAngle); 0],targetAnglesEMG(iAngle));
        
        cursorCir = circle(rCirCursor,EMGDatax,EMGDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if countBuffer == iterUpdatePlotEMG
            drawnow;
            countBuffer = 0;
        end
        
        if strcmp(state,tempState) && countState == 0
            countState = countState+1;
            xl = xlim;
            hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',24);
        elseif ~strcmp(state,tempState) && countState > 0
            countState = 0;
            delete(hsta)
        end

        tempState = state;
        
        % Trial
        switch state
            case 'start'
                xl = xlim; yl = ylim;
                delete(htrl)
                htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(trialNum)],'clipping','off','Fontsize',16);
                
                iAngle = iAngle+1;
                
                if reptarg
                    if iAngle == 1
                        iAngle = length(targetAnglesEMG);
                    else
                        iAngle = iAngle-1;
                    end
                    reptarg = 0;
                end
                if iAngle > length(targetAnglesEMG)
                    iAngle = 1;
                end
                
                %htrg = plot([0 targetPosx(iAngle)],[0 targetPosy(iAngle)],'r','Linewidth',3);
                targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                htrg = plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
                
                state = 'movement';
                tmove = tic;
              % When using undefined targets for EMGCO calibration
%             case 'movement'
%                 if toc(tmove) > movemtimeCal
%                     state = 'hold';
%                     tholdstart = tic;
%                 end
%             case 'hold'
%                 if toc(tholdstart) > holdtimeCal
%                     state = 'success';
%                     tsuccess = tic;
%                 end

      case 'movement'
                if cursorInTarget(cursorCir,targetCir) && toc(tmove) <= movemtimeCal
                    state = 'hold';
                    tholdstart = tic;
                elseif ~cursorInTarget(cursorCir,targetCir) && toc(tmove) > movemtimeCal
                    state = 'fail';
                    tfail = tic;
                end
            case 'hold'
                if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart) <= holdtimeCal
                    cursorHoldOut = cursorHoldOut+1;
                elseif cursorHoldOut >= 20  && toc(tholdstart) <= holdtimeCal
                    cursorHoldOut = 0;
                    state = 'fail';
                    tfail = tic;
                elseif cursorHoldOut > 20 && toc(tholdstart) > holdtimeCal
                    cursorHoldOut = 0;
                    state = 'fail';
                    tfail = tic;
                elseif cursorHoldOut <= 20 && toc(tholdstart) > holdtimeCal
                    cursorHoldOut = 0;
                    state = 'success';
                    tsuccess = tic;
                end
                
            case 'success'
                if toc(tsuccess) > timeout
                    state = 'relax';
                    trelax = tic;
                    delete(htrg)
                end
            case 'fail'
                if toc(tfail) > timeout
                    state = 'relax';
                    trelax = tic;
                    delete(htrg)
                    reptarg = 1;
                end
            case 'relax'
                if cursorInTarget(cursorCir,circle(2*rCirCursor,0,0)) && toc(trelax) > relaxtime
                    state = 'start';
                    trialNum = trialNum+1;
                end
        end
        
        % Appending trial data
        sampleNum = sampleNum+1;
        if saveforce
            forceDataOut_EMGCO(sampleNum,:) = {trialNum,iAngle,state,timeStamp,forceData(:,1),forceData(:,2),forceData(:,3),triggerData};
        end
        if saveEMG
            EMGDataOut_EMGCO(sampleNum,:) = {appendSamples'};
        end
    end
end