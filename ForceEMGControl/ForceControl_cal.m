%% DAQ from F/T Sensor and EMG - CO Task
% Virginia Casasnovas
% 11/07/2018

function ForceControl_cal(varargin)
%% Parameter assignment
% Default parameters
% File parameters
saveforce =     0;
saveEMG =       0;
date =          '20180711';
subject =       100;
task =          'CO';
code =          'calib';
filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
filepath =  [pwd '\Data\' date '\s' subject '\'];

% Task parameters
numTrials =         30;
targetForce =       5; % [N]
targetTol =         0.1;
cursorTol =         1.5;
numTargetsForce =   4;
targetAnglesForce = [pi/4:3*pi/(2*(numTargetsForce-1)):7*pi/4]; % [rad]
movemtime =         5; % sec
holdtime =          1; % sec
timeout =           1; % sec
relaxtime =         1; % sec

% Force parameters
fclF =              5; % [Hz]
scanRate =          2000; % [scans/sec]
availSamples =      40; % [samples]
bufferWin =         200; % [samples]
iterUpdatePlot =    10;
rotAngPlot =        0;

% EMG parameters
EMGEnabled = 0;
plotEMG =           0;
channelSubset =     [1 18];
channelName =       {'BB','Saw'};
sampleRateEMG =     1024;

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

%rCirTarget =        targetForce*targetTol; % [N]
rCirCursor =        targetForce*targetTol/cursorTol; % [N]

if strcmp(code,'calib')
    targetForce = 15;
    timeout = 0.5;
    movemtime = 2;
end

if length(channelSubset)~=length(channelName)
    error('Names for all channels not available.')
end

%% Initialization
disp('Running DataAcquisition for CO force task with EMG.')

device = daq.getDevices;

if ~isempty(device)
    disp('DAQ device found.')
    
    if saveEMG || saveforce
        if exist([filepath,filenameEMG],'file')
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
        

    
    % Initialize EMG
    if saveEMG
        disp('Initializing EMG.')
        library = TMSi.Library('usb');
        [EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channelSubset,channelName,sampleRateEMG);
    end
    
    if EMGEnabled
        disp('EMG initialized.')
        if plotEMG
            disp('Plotting EMGs to check.')
            EMGplot(emg_data,sampler,channels,channelSubset);
        end
    else
        disp('EMG could not be initialized.')
    end
    fprintf('\n')

        % Create NI DAQ session
    disp('Creating NI DAQ session.')
    s = daq.createSession('ni');
    addAnalogInputChannel(s,device.ID,0:6,'Voltage');
    if saveEMG
        addAnalogOutputChannel(s,device.ID,'ao0','Voltage');
    end
    
    % Obtain offset by averaging 2 sec of still data
    input('Press enter when prepared for sensor offset calculation.')
    fprintf('Obtaining offset values...\n')
    if saveEMG
        queueOutputData(s,zeros(2*scanRate,1));
    end
    forceOffsetData = s.startForeground();
    forceOffset = mean(forceOffsetData(:,1:6),1);
    disp('Offset obtained.')
    
    % Obtain Fz calibration by averaging 2 sec of holding data
    input('\nPress enter when holding handle.')
    fprintf('Obtaining Fz calibration values...\n')
    if saveEMG
        queueOutputData(s,zeros(2*scanRate,1));
    end
    calforceDataz = FzCalibration(s.startForeground,forceOffset);
    disp('Fz calibration obtained.')
    
    %% Data acquisition
    input('\nPress enter to start acquisition.')
    
    % Initialize variables
    global tmove trelax tfail tsuccess tholdstart
    global iAngle reptarg
    global htrg hsta htrl
    
    trialNum = 0;
    countState = 0;
    countBuffer = 0;
    forceDataBuffer = zeros(bufferWin,3);
    state = 'start';
    tempState = 'start';
    iAngle = 0;
    reptarg = 0;
    
    % Set target forces
    %targetAngles(targetAngles == 2*pi) = [];
    targetPosx = targetForce*cos(targetAnglesForce)';
    targetPosy = targetForce*sin(targetAnglesForce)';
    
    % Set figure
    hf = figure('Name','CO Force Control Task');
    [hf,hp] = Figinit(hf,[max(targetPosx) max(targetPosy)]./1.2);
    title('2D Force');
    xlabel('F_x [N]'); ylabel('F_y [N]');
    
    % Save file header
    if saveforce
        sampleNum = 1;
        forceDataOut_ForceCO(sampleNum,:) = {'Trialnum', 'TargetAng', 'State', 'TimeStamp', 'Fx', 'Fy', 'Fz','Trigger'};
    end
    if saveEMG
        sampleNum = 1;
        EMGDataOut_ForceCO(sampleNum,:) = {'EMG'};
    end
    % Start EMG data sampling
    if EMGEnabled
        sampler.start()
    end
    
    % Add event listeners and start acquisition
    if saveEMG
        outputData = [4*ones(1,availSamples),zeros(1,scanRate-availSamples)]';
        queueOutputData(s,outputData);
        hlout = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));
    end
    
    hlin = addlistener(s,'DataAvailable',@(src,event) processForceData(event,forceOffset,calforceDataz,hp));
    %hstop = addlistener(s,'DataAvailable',@(src,event) stopTrialNum(trialNum,numTrials));
    s.IsContinuous = true;
    s.Rate = scanRate; % scans/sec, samples/sec?
    s.NotifyWhenDataAvailableExceeds = availSamples; % Call listener when x samples are available
    s.startBackground();
    
    if trialNum == numTrials
        s.stop();
    else
        input('\Press enter to stop acquisition.')
    end
    
    % Save data
    if saveforce
        save([filepath,filenameforce],'forceDataOut_ForceCO')
    end
    if saveEMG && EMGEnabled
        save([filepath,filenameEMG],'EMGDataOut_ForceCO')
    end
     
    % Close session and delete handles
    s.stop()
    delete(hlin)
    if saveEMG
        delete(hlout)
    end
    %delete(hstop)
    close(hf)
    
    % Close EMG
    if saveEMG
        if EMGEnabled
            sampler.stop()
            sampler.disconnect()
        end
        library.destroy()
    end
    
else
    error('No DAQ device detected.')
end

%% Event handles
    function[forceDataz] = FzCalibration(event,offset)
        calibMat = [-0.56571    -0.01516    1.69417     -31.81016   -0.35339    33.73195;...
            -0.67022    37.27575    1.14848     -18.75980   0.34816     -19.33789;...
            19.06483    0.45530     18.57588    0.72627     19.51260    0.65624;...
            0.71302     0.12701     -32.12926   -1.43043    32.51061    1.09333;...
            37.80563    1.05190     -19.89516   -0.79823    -18.79634   -0.79836;...
            0.48287     -18.59045   0.38579     -18.05502   0.75146     -19.86235];
        offsetMat = repmat(offset,size(event(:,1:6),1),1);
        offsetData = event(:,1:6)-offsetMat;
        
        procData = zeros(size(event(:,1:6)));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        forceData = procData(:,1:3);      
        forceDataz = mean(forceData(:,3));
    end

    function processForceData(event,offset,calforceDataz,hp)
        calibMat = [-0.56571    -0.01516    1.69417     -31.81016   -0.35339    33.73195;...
            -0.67022    37.27575    1.14848     -18.75980   0.34816     -19.33789;...
            19.06483    0.45530     18.57588    0.72627     19.51260    0.65624;...
            0.71302     0.12701     -32.12926   -1.43043    32.51061    1.09333;...
            37.80563    1.05190     -19.89516   -0.79823    -18.79634   -0.79836;...
            0.48287     -18.59045   0.38579     -18.05502   0.75146     -19.86235];
        offsetMat = repmat(offset,size(event.Data(:,1:6),1),1);
        offsetData = event.Data(:,1:6)-offsetMat;
        
        procData = zeros(size(event.Data(:,1:6)));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        timeStamp = event.TimeStamps;
        triggerData = event.Data(:,7);
        forceData = procData(:,1:3);
        
        bufferTemp = forceDataBuffer;
        bufferTemp(1:availSamples,:) = forceData;
        bufferTemp(availSamples+1:end,:) = forceDataBuffer(1:bufferWin-availSamples,:);
        forceDataBuffer = bufferTemp;
        countBuffer = countBuffer+1;
        
        % LPF
        wn = (2/s.Rate)*fclF;
        [b,a] = butter(2,wn,'low');
        forceDataFilt = filter(b,a,forceDataBuffer);
        
        forceDatax = mean(forceDataFilt(:,1));
        forceDatay = mean(forceDataFilt(:,2));
        forceDataz = mean(forceDataFilt(:,3));
        
        calCirCursor = rCirCursor*(1 + abs(forceDataz/calforceDataz));
        forceDataRot = Rot([forceDatax,forceDatay], rotAngPlot);
        cursorCir = circle(calCirCursor,forceDataRot(1),forceDataRot(2));
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if countBuffer == iterUpdatePlot
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
               
        % Trial
        tempState = state;
        
        switch state
            case 'start'
                xl = xlim; yl = ylim;
                delete(htrl)
                htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(trialNum)],'clipping','off','Fontsize',16);
       
                iAngle = iAngle+1;
                
                if reptarg
                    if iAngle == 1
                        iAngle = length(targetAnglesForce);
                    else
                        iAngle = iAngle-1;
                    end
                    reptarg = 0;
                end
                if iAngle > length(targetAnglesForce)
                    iAngle = 1;
                end
                
                htrg = plot([0 targetPosx(iAngle)],[0 targetPosy(iAngle)],'r','Linewidth',3);
                
                
                state = 'movement';
                tmove = tic;
            case 'movement'
                if calCirCursor > 2*rCirCursor
                    state = 'fail';
                    tfail = tic;
                else
                    if toc(tmove) > movemtime
                        state = 'hold';
                        tholdstart = tic;
                    end
                end
            case 'hold'
                if calCirCursor > 2*rCirCursor
                    state = 'fail';
                    tfail = tic;
                else
                    if toc(tholdstart) > holdtime
                        state = 'success';
                        tsuccess = tic;
                    end
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
                if cursorInTarget(cursorCir,circle(1.5*calCirCursor,0,0)) && toc(trelax) > relaxtime && calCirCursor <= 2*rCirCursor
                    state = 'start';
                    trialNum = trialNum+1;
                end
        end
        
        % Appending trial data 
        if saveforce 
            sampleNum = sampleNum+1;
            forceDataOut_ForceCO(sampleNum,:) = {trialNum,iAngle,state,timeStamp,forceData(:,1),forceData(:,2),forceData(:,3),triggerData};
        end
        % Appending EMG data
        if saveEMG && EMGEnabled
            samples = sampler.sample();
            appendSamples = samples(channelSubset,:);
            EMGDataOut_ForceCO(sampleNum,:) = {appendSamples'};
        end
    end

%     function stopTrialNum(trialNum,numTrials)
%         if trialNum == numTrials
%             s.stop();
%             close(hf)
%         end
%     end
end
