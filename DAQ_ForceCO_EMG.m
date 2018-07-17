%% DAQ from F/T Sensor and EMG - CO Task
% Virginia Casasnovas
% 11/07/2018

function DAQ_ForceCO_EMG(varargin)
%% Parameter assignment
% Default parameters
% File parameters
saveforce =  0;
saveEMG =    0;
date =      '20180711';
task =      'CO';
code =      '001';
filenameforce =  [date,'_',task,'_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filepath =  [pwd '\Data\' date '\'];

% Task parameters
targetForce =   5; % [N]
numTargets =    8;
rCirTarget =    targetForce/10; % [N]
rCirCursor =    targetForce/20; % [N]

scanRate =      1000; % [scans/sec]
availSamples =  100; % [samples]
fc =            5; % [Hz]

movemtime =     5; % sec
holdtime =      1; % sec
timeout =       1; % sec
relaxtime =     1; % sec

% EMG parameters
plotEMG =           0;
EMGEnabled =        0;
channel_subset =    [1 18];
channel_name =      {'BB','Saw'};
sample_rate =       1024;

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

if length(channel_subset)~=length(channel_name)
    error('Names for all channels not available.')
end

%% Initialization
disp('Running DataAcquisition for CO force task with EMG.')
device = daq.getDevices;

if ~isempty(device)
    if saveforce
        if exist([filepath,filenameforce],'file')
            savefile = input('\nThis file already exsists. Continue saving? (y/n) ','s');
            if strcmp(savefile,'y')
                fprintf('Saving force data in %s.\n',filenameforce)
                if saveEMG
                    fprintf('Saving EMG data in %s.\n\n',filenameEMG)
                else
                    fprintf('\n')
                end
            else
                saveforce = 0;
                fprintf('Not saving data.\n\n')
            end
        else
            fprintf('Saving force data in %s.\n',filenameforce)
            if saveEMG
                fprintf('Saving EMG data in %s.\n\n',filenameEMG)
            else
                fprintf('\n')
            end
        end    
    else
        fprintf('Not saving data.\n\n')
    end
        
    % Create NI DAQ session
    disp('Creating NI DAQ session.')
    s = daq.createSession('ni');
    addAnalogInputChannel(s,device.ID,0:6,'Voltage');
    addAnalogOutputChannel(s,device.ID,'ao0','Voltage');
    
    % Initialize EMG
    disp('Initializing EMG.')
    library = TMSi.Library('usb');
    [EMGEnabled,sampler,emg_data,channels] = EMGinit(library,channel_subset,channel_name,sample_rate);
    
    if EMGEnabled
        disp('EMG initialized.')
        if plotEMG
            disp('Plotting EMGs to check.')
            emg_plot = TMSi.RealTimePlot('EMG RealTimePlot', sampler.sample_rate, channels);
            emg_plot.setWindowSize(10);
            emg_plot.show();
            
            sampler.start();
            samples = sampler.sample();
            emg_data.append(samples(channel_subset,:));
            while emg_plot.is_visible
                emg_plot.append(samples(channel_subset,:));
                emg_plot.draw();
            end
            sampler.stop();
        end
    else
        disp('EMG could not be initialized.')
    end
    
    % Obtain offset by averaging 2 sec of still data
    fprintf('\nObtaining offset values...\n')
    queueOutputData(s,zeros(2*scanRate,1));
    forceOffsetData = s.startForeground();
    forceOffset = mean(forceOffsetData(:,1:6),1);
    disp('Offset obtained.')
    
    % Obtain Fz calibration by averaging 2 sec of holding data
    input('\nPress enter when holding knob.')
    fprintf('Obtaining Fz calibration values...\n')
    queueOutputData(s,zeros(2*scanRate,1));
    calforceDataz = FzCalibration(s.startForeground,forceOffset);
    disp('Fz calibration obtained.')
    
    % Start data acquisition
    input('\nPress enter to start acquisition.')
    
    % Set target forces
    targetAngles = [0:2*pi/numTargets:2*pi]; % [rad]
    targetAngles(targetAngles == 2*pi) = [];
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    
    % Set figure
    hf = figure('Name','CO Force Control Task');
    [hf,hp] = Figinit(hf,target);
    title('2D Force');
    xlabel('F_x [N]'); ylabel('F_y [N]');
    
    % Set file if saving option is enabled
    if saveforce
        samplenum = 1;
        forceDataOut(samplenum,:) = {'Trialnum', 'State', 'TimeStamp', 'Fx', 'Fy', 'Trigger'};
    end
    
    %% Data acquisition
    % Initialize variables
    global tstart trelax tfail tsuccess tholdstart
    global targetCir
    global htrg hsta htrl
    
    trialnum = 0;
    cursorholdout = 0;
    countstate = 0;
    state = 'start';
    prevstate = 'start';
    
    % Start EMG data sampling
    if EMGEnabled
        sampler.start()
    end
    
    % Add event listener and start acquisition
    outputData = [4*ones(1,availSamples),zeros(1,scanRate-availSamples)]';
    queueOutputData(s,outputData);
    hlout = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));

    hlin = addlistener(s,'DataAvailable',@(src,event) processForceData(event,forceOffset,calforceDataz,hp));
    s.IsContinuous = true;
    s.Rate = scanRate; % scans/sec, samples/sec?
    s.NotifyWhenDataAvailableExceeds = availSamples; % Call listener when x samples are available
    s.startBackground();
    
    input('\Press enter to stop acquisition.') 
    
    % Save data
    if saveforce
        save([filepath,filenameforce],'forceDataOut')
    end
    if saveEMG && EMGEnabled
        EMGDataOut = emg_data.samples;
        save([filepath,filenameEMG],'EMGDataOut')
    end
     
    % Close session and delete handles
    s.stop()
    delete(hlin)
    delete(hlout)
    close(hf)
    
    % Close EMG
    if EMGEnabled
        sampler.stop()
        sampler.disconnect()
    end
    library.destroy()
    
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
        
%         % LPF?
%         wn = (2/s.Rate)*fc;
%         [b,a] = butter(2,wn,'low');
%         forceData = filter(b,a,forceData);
        
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
        
        % LPF?
        wn = (2/s.Rate)*fc;
        [b,a] = butter(2,wn,'low');
        forceDataFilt = filter(b,a,forceData);
        
        fmx = mean(forceData(:,1));
        fmy = mean(forceData(:,2));
        forceDatax = mean(forceDataFilt(:,1));
        forceDatay = mean(forceDataFilt(:,2));
        forceDataz = mean(forceDataFilt(:,3));
        fprintf('[%f,%f]   [%f,%f]\n',fmx,fmy,forceDatax,forceDatay)
        calCirCursor = rCirCursor*(1 + abs(forceDataz/calforceDataz));
        cursorCir = circle(calCirCursor,forceDatax,forceDatay);
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
                if calCirCursor > rCirCursor*2
                    state = 'fail';
                    tfail = tic;
                else
                    if cursorInTarget(cursorCir,targetCir) && toc(tstart) <= movemtime
                        state = 'hold';
                        tholdstart = tic;
                    elseif ~cursorInTarget(cursorCir,targetCir) && toc(tstart) > movemtime
                        state = 'fail';
                        tfail = tic;
                    end
                end
            case 'hold'
                if calCirCursor > rCirCursor*2
                    state = 'fail';
                    tfail = tic;
                else
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
                if cursorInTarget(cursorCir,circle(1.5*calCirCursor,0,0)) && toc(trelax) > relaxtime && calCirCursor <= rCirCursor*2
                    state = 'start';
                    trialnum = trialnum+1;
                end
        end
        
        % Saving EMG data in emg_data matrix
        if EMGEnabled
            samples = sampler.sample();
            emg_data.append(samples(channel_subset,:));
        end
        
        % Saving trial data in cell
        if saveforce 
            samplenum = samplenum+1;
            forceDataOut(samplenum,:) = {trialnum,state,timeStamp,forceData(:,1),forceData(:,2),triggerData};
        end
    end
end
