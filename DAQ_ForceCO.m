%% DAQ from F/T Sensor - CO Task
% Virginia Casasnovas
% 07/07/2018

function DAQ_ForceCO(varargin)
%% Parameter assignment

% Default parameters
% File parameters
saveforce =  0;
date =      '20180709';
task =      'CO';
code =      '001';
filenameforce =  [date,'_',task,'_',code,'.txt'];
filepath =  [pwd '\Data\' date '\'];

% Task parameters
targetForce =   5; % [N]
numTargets =    8;
rCirTarget =    targetForce/10; % [N]
rCirCursor =    targetForce/20; % [N]

scanRate =      1000; % [scans/sec]
availSamples =  100; % [samples]
fc =            8; % [Hz]

movemtime =     5; % sec
holdtime =      1; % sec
timeout =       1; % sec
relaxtime =     1; % sec

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(who,varargin{ii});
    end
end

%% Initialization
disp('Running DataAcquisition for CO force task.')
device = daq.getDevices;


if ~isempty(device)
    if saveforce
        if exist([filepath,filenameforce],'file')
            savefile = input('This file already exsists. Continue saving? (y/n) ','s');
            if strcmp(savefile,'y')
                fprintf('Saving data in %s.\n\n',filenameforce)
            else
                saveforce = 0;
                disp('Not saving data in txt file.')
            end
        end
    else
        disp('Not saving data in txt file.')
    end
    
    disp('Creating session.')
    s = daq.createSession('ni');
    addAnalogInputChannel(s,device.ID,0:5,'Voltage');
    
    % Obtain offset by averaging 2 sec of still data
    fprintf('\nObtaining offset values...\n')
    s.DurationInSeconds = 2;
    offset = mean(s.startForeground(),1);
    disp('Offset obtained.')
    
    % Obtain Fz calibration by averaging 2 sec of holding data
    input('\nPress enter when holding knob.')
    fprintf('Obtaining Fz calibration values...\n')
    s.DurationInSeconds = 2;
    calforceDataz = FzCalibration(s.startForeground,offset);
    disp('Fz calibration obtained.')
    
    % Set target forces
    targetAngles = [0:2*pi/numTargets:2*pi]; % [rad]
    targetAngles(targetAngles == 2*pi) = [];
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    
    % Start data acquisition
    input('\nPress enter to start acquisition.')
    
    % Set figure
    hf = figure('Name','CO Force Control Task');
    set(hf,'units','normalized','outerposition',[0 0 1 1]);
    set(hf,'menubar','none')
    set(hf,'color','w');
    hp = patch(0,0,'g','Edgecolor','g');
    hold on;
    title('2D Force');
    xlabel('F_x [N]'); ylabel('F_y [N]');
    set(gca,'xlim',targetForce*[-1.2 1.2],'ylim',targetForce*[-1.2 1.2]);
    set(gca,'TickLength',[0 0])
    axis square;
    grid on;
    box on;
    
    % Set file if saving option is enabled
    if saveforce
        fid = fopen([filepath,filenameforce],'wt');
        fprintf(fid, 'Trialnum\t State\t TimeStamp\t Fx [N]\t Fy [N]\n');
    end
    
    %% Data acquisition
    % Initialize variables
    global tstart tmovem trelax tfail tsuccess tholdstart
    global targetCir
    global htrg hsta htrl
    
    trialnum = 0;
    cursorholdout = 0;
    countstate = 0;
    state = 'start';
    prevstate = 'start';
    
    % Add event listener and start acquisition
    hl = addlistener(s,'DataAvailable',@(src,event) processData(event,offset,calforceDataz,hp));
    s.IsContinuous = true;
    s.Rate = scanRate; % scans/sec, samples/sec?
    s.NotifyWhenDataAvailableExceeds = availSamples; % Call listener when x samples are available
    s.startBackground();
    
    input('\Press enter to stop acquisition.')
    
    % Close session and delete handles
    s.stop()
    delete(hl)
    close(hf)
    if saveforce
        fclose(fid);
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
        offsetMat = repmat(offset,size(event,1),1);
        offsetData = event-offsetMat;
        procData = zeros(size(event));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        forceData = procData(:,1:3);
        
        % LPF?
        %fc = 8;
        wn = (2/s.Rate)*fc;
        [b,a] = butter(3,wn,'low');
        forceData = filter(b,a,forceData);
        
        forceDataz = mean(forceData(:,3));
    end

    function processData(event,offset,calforceDataz,hp)
        calibMat = [-0.56571    -0.01516    1.69417     -31.81016   -0.35339    33.73195;...
            -0.67022    37.27575    1.14848     -18.75980   0.34816     -19.33789;...
            19.06483    0.45530     18.57588    0.72627     19.51260    0.65624;...
            0.71302     0.12701     -32.12926   -1.43043    32.51061    1.09333;...
            37.80563    1.05190     -19.89516   -0.79823    -18.79634   -0.79836;...
            0.48287     -18.59045   0.38579     -18.05502   0.75146     -19.86235];
        offsetMat = repmat(offset,size(event.Data,1),1);
        offsetData = event.Data-offsetMat;
        procData = zeros(size(event.Data));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        timeStamp = event.TimeStamps;
        forceData = procData(:,1:3);
        
        % LPF?
        %fc = 8;
        wn = (2/s.Rate)*fc;
        [b,a] = butter(2,wn,'low');
        forceData = filter(b,a,forceData);
        
        forceDatax = mean(forceData(:,1));
        forceDatay = mean(forceData(:,2));
        forceDataz = mean(forceData(:,3));
        calCirCursor = rCirCursor*(1 + abs(forceDataz/calforceDataz));
        cursorCir = circle(calCirCursor,forceDatax,forceDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        drawnow;

        if strcmp(state,prevstate) && countstate == 0
            countstate = countstate+1;
            xl = xlim;
            hsta = text(xl(2)+0.3*xl(2),0,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
        elseif ~strcmp(state,prevstate) && countstate > 0
            countstate = 0;
            delete(hsta)
        end
        
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
                        tmovem = toc(tstart);
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
        
        % Saving data in txt file
        if saveforce
            fprintf(fid, '%d\t %s\t %1.4f\t %1.4f\t %1.4f\n',trialnum,state,timeStamp(1),forceDatax,forceDatay);
        end
    end
end
