%% DAQ from EMG - CO Task
% Virginia Casasnovas
% 11/07/2018

function DAQ_EMGCO(varargin)
%% Parameter assignment
% Default parameters
% File parameters
saveforce =  0;
saveEMG =    1;
device =     [];
date =      '20180717';
task =      'CO';
code =      '001';
filenameforce = [date,'_',task,'_F_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filepath =  [pwd '\Data\' date '\'];

% Task parameters
num_trials =    3;
targetForce =   1000; % [N]
numTargets =    3;
rCirTarget =    targetForce/10; % [N]
rCirCursor =    targetForce/20; % [N]

fc =            10; % [Hz]

scanRate =      1000; % [scans/sec]
availSamples =  100; % [samples]

movemtime =     5; % sec
holdtime =      1; % sec
timeout =       1; % sec
relaxtime =     1; % sec

% EMG parameters
plotEMG =           0;
EMGEnabled =        0;
channel_subset =    [1 2];
channel_name =      {'BB','TL'};
sample_rate =       1024;
smooth_win =        200;

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

if length(channel_subset)~=length(channel_name)
    error('Names for all channels not available.')
end

disp('Running DataAcquisition for CO EMG task.')

if saveEMG
    if exist([filepath,filenameEMG],'file')
        savefile = input('\nThis file already exsists. Continue saving? (y/n) ','s');
        if strcmp(savefile,'y')
            fprintf('Saving EMG data in %s.\n',filenameEMG)
            if saveEMG
                fprintf('Saving force data in %s.\n\n',filenameforce)
            else
                fprintf('\n')
            end
        else
            saveEMG = 0;
            fprintf('Not saving data.\n\n')
        end
    else
        fprintf('Saving force data in %s.\n',filenameEMG)
        if saveforce
            fprintf('Saving EMG data in %s.\n\n',filenameforce)
        else
            fprintf('\n')
        end
    end
else
    fprintf('Not saving data.\n\n')
end

% Save header
if saveEMG
    samplenum = 1;
    EMGDataOut(samplenum,:) = {'Trialnum', 'State', 'EMG'};
end

%% Initialize EMG
disp('Initializing EMG.')
library = TMSi.Library('usb');
[EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channel_subset,channel_name,sample_rate);

if EMGEnabled
    disp('EMG initialized.')
    
    if plotEMG
        disp('Plotting EMGs to check.')
        EMGplot(emg_data,sampler,channels,channel_subset);
    end
    
    %% Initialize force DAQ
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
            
            % Header for forcefile
            if saveforce
                samplenum = 1;
                forceDataOut(samplenum,:) = {'Trialnum', 'State', 'TimeStamp', 'Force', 'Trigger'};
            end
        else
            device = [];
            disp('DAQ device not found.')
        end
    end
    %% Start data acquisition
    input('\nPress enter to start acquisition.')
    
    % Set target forces
    targetAngles = [pi/4:pi/(2*(numTargets-1)):3*pi/4]; % [rad]
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    
    % Set figure
    global htrl
    hf = figure('Name','CO EMG Control Task');
    [hf,hp] = Figinit(hf,targetForce);
    title('EMG');
    xlabel(channel_name{1}); ylabel(channel_name{2});
    xl = xlim; yl = ylim;
    htrl = text(xl(2)+0.3*xl(2),yl(2),'Trial: 0','clipping','off','Fontsize',14);
    %set(gca,'XTickLabel',[],'YTickLabel',[])
    
    %% Data acquisition
    % Start EMG data sampling
    sampler.start()
    countstate = 0;
    
    if ~isempty(device)
        % Initialize variables
        global tstart trelax tfail tsuccess tholdstart
        global targetCir
        global htrg hsta
        
        trialnum = 0;
        cursorholdout = 0;
        countstate = 0;
        state = 'start';
        prevstate = 'start';
        
        % Add event listener and start acquisition
        outputData = [4*ones(1,availSamples),zeros(1,scanRate-availSamples)]';
        queueOutputData(s,outputData);
        hlout = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));
        
        hlin = addlistener(s,'DataAvailable',@(src,event) processForceData(event,forceOffset,hp));
        s.IsContinuous = true;
        s.Rate = scanRate; % scans/sec, samples/sec?
        s.NotifyWhenDataAvailableExceeds = availSamples; % Call listener when x samples are available
        s.startBackground();

        input('\Press enter to stop acquisition.') 

    else
        disp('Device not found');
        for itrial = 1:num_trials
            
            nexttrial = false;
            cursorholdout = 0;
            state = 'start';
            prevstate = 'start';
            
            while ~nexttrial
                samples = sampler.sample();
                samples_subset = samples(channel_subset,:);
                tvis = tic;
                pause(0.02)
                if toc(tvis) >= 0.02
                    emg_data.append(samples_subset);
                    emg_data_len = size(emg_data.samples,2);
                    % fprintf('%d\n',size(samples_subset,2));
                    % LPF?
                    wn = (2/sample_rate)*fc;
                    [b,a] = butter(2,wn,'high');
                    filtemg = filter(b,a,emg_data.samples);
                    emg_proc = mean(abs(filtemg),2); % Rectify and average
                    
                    [EMGDatax,EMGDatay] = EMG2xy(emg_proc);
                    
                    cursorCir = circle(rCirCursor,EMGDatax,EMGDatay);
                    set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
                    
                    if strcmp(state,prevstate) && countstate == 0
                        countstate = countstate+1;
                        xl = xlim;
                        hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
                    elseif ~strcmp(state,prevstate) && countstate > 0
                        countstate = 0;
                        delete(hsta)
                    end
                    
                    drawnow;
                    
                    prevstate = state;
                    
                    % Trial
                    switch state
                        case 'start'
                            tstart = tic;
                            
                            xl = xlim; yl = ylim;
                            delete(htrl)
                            htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(itrial)],'clipping','off','Fontsize',14);
                            
                            iAngle = randi(numTargets);
                            targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                            htrg = plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
                            
                            state = 'movement';
                        case 'movement'
                            if cursorInTarget(cursorCir,targetCir) && toc(tstart) <= movemtime
                                state = 'hold';
                                tholdstart = tic;
                            elseif ~cursorInTarget(cursorCir,targetCir) && toc(tstart) > movemtime
                                state = 'fail';
                                tfail = tic;
                            end
                        case 'hold'
                            if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart) <= holdtime
                                cursorholdout = cursorholdout+1;
                            elseif cursorholdout > 5 && toc(tholdstart) > holdtime
                                cursorholdout = 0;
                                state = 'fail';
                                tfail = tic;
                            elseif cursorholdout <= 5 && toc(tholdstart) > holdtime
                                cursorholdout = 0;
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
                                nexttrial = true;
                            end
                    end
                end
            end
        end
        
        if saveEMG
            samplenum = samplenum+1;
            EMGDataOut(samplenum,:) = {itrial,state,samples_subset'};
        end
    end
     
    % Close EMG
    sampler.stop()
    sampler.disconnect()
    
    % Delete handles and stop session
    close(hf)
    if ~isempty(device)
        s.stop()
        delete(hlin)
        delete(hlout)
    end
    
    % Save mat
    if saveEMG
        save([filepath,filenameEMG],'EMGDataOut')
    end
else
    disp('EMG could not be initialized.')
end
library.destroy()
%% Function handles
    function processForceData(event,offset,hp)
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
        
        samples = sampler.sample();
        samples_subset = samples(channel_subset,:);
        emg_data.append(samples_subset);
        emg_data_len = size(emg_data.samples,2);

        wn = (2/sample_rate)*fc;
        [b,a] = butter(2,wn,'high');
        filtemg = filter(b,a,emg_data.samples);
        emg_proc = mean(abs(filtemg),2); % Rectify and average
        
        [EMGDatax,EMGDatay] = EMG2xy(emg_proc);
        
        cursorCir = circle(rCirCursor,EMGDatax,EMGDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if strcmp(state,prevstate) && countstate == 0
            countstate = countstate+1;
            xl = xlim;
            hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
        elseif ~strcmp(state,prevstate) && countstate > 0
            countstate = 0;
            delete(hsta)
        end
        
        drawnow;
        
        prevstate = state;
        
        % Trial
        switch state
            case 'start'
                tstart = tic;
                
                xl = xlim; yl = ylim;
                delete(htrl)
                htrl = text(xl(2)+0.3*xl(2),yl(2),['Trial: ',num2str(trialnum)],'clipping','off','Fontsize',14);
                
                iAngle = randi(numTargets);
                targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                htrg = plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
                
                state = 'movement';
            case 'movement'
                if cursorInTarget(cursorCir,targetCir) && toc(tstart) <= movemtime
                    state = 'hold';
                    tholdstart = tic;
                elseif ~cursorInTarget(cursorCir,targetCir) && toc(tstart) > movemtime
                    state = 'fail';
                    tfail = tic;
                end
            case 'hold'
                if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart) <= holdtime
                    cursorholdout = cursorholdout+1;
                elseif cursorholdout > 5 && toc(tholdstart) > holdtime
                    cursorholdout = 0;
                    state = 'fail';
                    tfail = tic;
                elseif cursorholdout <= 5 && toc(tholdstart) > holdtime
                    cursorholdout = 0;
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
                    trialnum = trialnum + 1;
                end
        end
                
        % Saving trial data in cell
        if saveforce
            samplenum = samplenum+1;
            forceDataOut(samplenum,:) = {trialnum,state,timeStamp,forceData,triggerData};
        end
    end
end