%% DAQ from F/T Sensor - Bar Task
% Virginia Casasnovas
% 07/07/2018

function DAQ_ForceBAR(varargin)
%% Parameter assignment

% Default parameters
% File parameters
saveforce =  0;
date =      '20180709';
task =      'BAR';
code =      '001';
filename =  [date,'_',task,'_',code,'.txt'];
filepath =  [pwd '\Data\',filename];

% Task parameters
targetForce =   10; % [N]
barTol =        targetForce/10; % [N]

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
disp('Running DataAquisition for bar force task.')
device = daq.getDevices;

if ~isempty(device)
    if saveforce
        if exist(filepath,'file')
            savefile = input('This file already exsists. Continue saving? (y/n) ','s');
            if strcmp(savefile,'y')
                fprintf('Saving data in %s.\n',filename)
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
    
    % Set target forces
    targetBar = targetForce;
    
    % Start data acquisition
    input('\nPress enter to start acquisition.')
    
    % Set figure
    hf = figure('Name','Bar Force Control Task');
    set(hf,'units','normalized','outerposition',[0 0 1 1]);
    set(hf,'menubar','none')
    set(hf,'color','w');
    bar(targetForce,'FaceAlpha',0.1,'BarWidth',0.3)
    hold on;
    hp = bar(0,'BarWidth',0.3);
    title('Force');
    ylabel('F [N]'); %xlabel('Target');
    set(gca,'ylim',targetForce*[0 1.2]);
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',[]);
    box on;
    axis square;
    
    % Set file if saving option is enabled
    if saveforce
        fid = fopen(filepath,'wt');
        fprintf(fid, 'Trialnum\t State\t TimeStamp\t Fx [N]\t Fy [N]\n');
    end
    
    %% Data acquisition
    % Initialize variables
    global tstart tmovem trelax tfail tsuccess tholdstart
    global htxt
    
    trialnum = 0;
    cursorholdout = 0;
    countstate = 0;
    state = 'start';
    prevstate = 'start';
    
    % Add event listener and start acquisition
    h = addlistener(s,'DataAvailable',@(src,event) processData(event,offset,hp));
    s.IsContinuous = true;
    s.Rate = 1000; % scans/sec, samples/sec?
    %s.NotifyWhenDataAvailableExceeds = 100; % Call listener when x samples are available
    s.startBackground();
    
    input('\Press enter to stop acquisition.')
    
    % Close session and delete handles
    s.stop()
    delete(h)
    close(hf)
    if saveforce
        fclose(fid);
    end
else
    error('No DAQ device detected.')
end

%% Event handles
    function processData(event,offset,hp)
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
        forceData = procData(:,1:2);
        
        %         % LPF?
        %         fc = 250;
        %         wn = (2/s.Rate)*fc;
        %         [b,a] = butter(2,wn,'low');
        %         forceData = filter(b,a,forceData);
        
        forceDatax = mean(forceData(:,1));
        forceDatay = mean(forceData(:,2));
        cursorBar = sqrt(forceDatax^2+forceDatay^2);
        set(hp,'ydata',cursorBar);
        
        if strcmp(state,prevstate) && countstate == 0
            countstate = countstate+1;
            xl = xlim; yl = ylim;
            htxt = text(xl(2)+0.1*xl(2),yl(2)/2,[upper(state(1)),state(2:end)],'clipping','off','Fontsize',16);
        elseif ~strcmp(state,prevstate) && countstate > 0
            countstate = 0;
            delete(htxt)
        end
        
        drawnow;
        
        prevstate = state;
        
        % Trial
        switch state
            case 'start'
                trialnum = trialnum+1;
                tstart = tic;
                state = 'movement';
            case 'movement'
                if (cursorBar <= targetBar+barTol && cursorBar >= targetBar-barTol) && toc(tstart) <= movemtime
                    tmovem = toc(tstart);
                    state = 'hold';
                    tholdstart = tic;
                elseif ~(cursorBar <= targetBar+barTol && cursorBar >= targetBar-barTol) && toc(tstart) > movemtime
                    state = 'fail';
                    tfail = tic;
                end
            case 'hold'
                if ~(cursorBar <= targetBar+barTol && cursorBar >= targetBar-barTol) && toc(tholdstart) <= holdtime
                    cursorholdout = cursorholdout+1;
                elseif cursorholdout > 10 && toc(tholdstart) > holdtime
                    cursorholdout = 0;
                    state = 'fail';
                    tfail = tic;
                elseif cursorholdout <= 10 && toc(tholdstart) > holdtime
                    cursorholdout = 0;
                    state = 'success';
                    tsuccess = tic;
                end
            case 'fail'
                if toc(tfail) > timeout
                    state = 'relax';
                    trelax = tic;
                end
            case 'success'
                if toc(tsuccess) > timeout
                    state = 'relax';
                    trelax = tic;
                end
            case 'relax'
                if cursorBar <= targetBar/10 && toc(trelax) > relaxtime
                    state = 'start';
                end
        end
        
        % Saving data in txt file
        if saveforce
            fprintf(fid, '%d\t %s\t %1.4f\t %1.4f\t %1.4f\n',trialnum,state,timeStamp(1),forceDatax,forceDatay);
        end
    end
end