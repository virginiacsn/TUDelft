%% DAQ from F/T Sensor
% Virginia Casasnovas
% 07/07/2018

function DataAcquisition(varargin)
%% Parameter assignment

% Default parameters
% File parameters
savefile = 0;
filename = '20180706_FS_test.txt';
filepath = ['D:\Student_experiments\Virginia\FTSensor\Data\',filename];

% Task parameters
targetForce = 8; % [N]
numTargets = 8; 
targetAngles = [0:2*pi/numTargets:2*pi]; % [rad]
targetAngles(targetAngles == 2*pi) = [];
rCirTarget = 0.5; % [N]
rCirCursor = 0.25; % [N]

movemtime = 5; % sec
holdtime = 0.5; % sec
timeout = 1; % sec
relaxtime = 1; % sec

% Overwrite parameters from input structs
if ~isempty(varargin)
    for ii = 1:length(varargin)
        struct2vars(varargin{ii});
    end
end

%% Initialization

device = daq.getDevices;

if ~isempty(device)
    s = daq.createSession('ni');
    addAnalogInputChannel(s,device.ID,0:5,'Voltage');
    
    % Obtain offset by averaging 1sec of still data
    disp('Obtaining offset values...')
    offset = mean(s.startForeground(),1);
    disp('Offset obtained.')
    
    input('\nPress enter to start acquisition.')
    
    % Set target forces
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    
    % Set figure
    hf = figure(1);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %scrsize = get(0, 'Screensize');
    %set(gcf, 'Position', [scrsize(3:4)/4 scrsize(3:4)/2]);
    set(hf,'units','normalized','outerposition',[0 0 1 1]);
    set(hf,'menubar','none')
    set(hf,'color','w');
    hp = plot(0,0,'g','LineWidth',2);
    hold on;
    %plot(targetPosx,targetPosy,'r*')
    %     for j = 1:length(targetPosx)
    %         targetCir = circle(rCirTarget,targetPosx(j),targetPosy(j));
    %         plot(targetCir(:,1),targetCir(:,2),'r','LineWidth',2)
    %         %rectangle('Position',[targetPosx(j),targetPosy(j),rCirTarget,rCirTarget],...
    %         %    'Curvature',[1,1],'FaceColor','r','EdgeColor','r');
    %     end
    title('2D Force');
    xlabel('F_x [N]'); ylabel('F_y [N]');
    set(gca,'xlim',targetForce*[-1.2 1.2],'ylim',targetForce*[-1.2 1.2]);
    axis square;
    grid on;
    
    % Set file if saving option is enabled
    if savefile
        fid = fopen(filepath,'wt');
        fprintf(fid, 'Trialnum\t State\t TimeStamp\t Fx [N]\t Fy [N]\n');
    end
    
    %% Data acquisition
    
    trialnum = 0;
    cursorholdout = 0;
    countstate = 0;
    state = 'start';
    prevstate = 'start';
    
    global tstart tmovem trelax tfail tsuccess tholdstart
    global targetCir
    global ht htxt
    
    h = addlistener(s,'DataAvailable',@(src,event) processData(event,offset,hp));
    s.IsContinuous = true;
    s.Rate = 1000; % scans/sec, samples/sec?
    %s.NotifyWhenDataAvailableExceeds = 100; % Call listener when x samples are available
    s.startBackground();
    
    input('\Press enter to stop acquisition.')
    
    s.stop()
    delete(h)
    close(hf)
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
        forceDatax = mean(forceData(:,1));
        forceDatay = mean(forceData(:,2));
        cursorCir = circle(rCirCursor,forceDatax,forceDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if strcmp(state,prevstate) && countstate == 0
            countstate = countstate+1;
            fprintf('%d ',countstate)
            htxt = text(6.5, 4, state, 'clipping', 'off','Fontsize',12);
        elseif ~strcmp(state,prevstate) && countstate > 0
            countstate = 0;
            delete(htxt)
        end
        
        drawnow;
        
        prevstate = state;
        
        % Trial epochs
        switch state
            case 'start'
                trialnum = trialnum+1;
                tstart = tic;
                
                iAngle = randi(8);
                targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
                ht = plot(targetCir(:,1),targetCir(:,2),'r','LineWidth',2);
                
                state = 'movement';
            case 'movement'
                if cursorInTarget(cursorCir,targetCir) && (toc(tstart)-tstart) <= movemtime
                    tmovem = toc(tstart);
                    state = 'hold';
                    tholdstart = tic;
                elseif (toc(tstart)-tstart)> movemtime
                    state = 'fail';
                end
            case 'hold'
                if ~cursorInTarget(cursorCir,targetCir) && toc(tholdstart)-tmovem <= holdtime
                    cursorholdout = cursorholdout+1;
                elseif cursorholdout > 10 && toc(tholdstart)-tmovem > holdtime
                    state = 'fail';
                    tfail = tic;
                elseif cursorholdout <= 10 && toc(tholdstart)-tmovem > holdtime
                    state = 'success';
                    tsuccess = tic;
                end
            case 'fail'
                if toc(tfail) > timeout
                    state = 'relax';
                    delete(ht)
                end
            case 'success'
                if toc(tsuccess) > timeout
                    state = 'relax';
                    delete(ht)
                end
            case 'relax'
                if cursorInTarget(cursorCir,circle(rCirCursor+0.1,0,0)) && toc(trelax) > relaxtime
                    state = 'start';
                end
        end
        
        if savefile
            fprintf(fid, '%d\t %s\t %1.4f\t %1.4f\t %1.4f\n',trialnum,state,timeStamp(1),forceDatax,forceDatay);
        end
    end
end
