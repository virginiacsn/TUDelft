% DAQ from F/T Sensor
function DataAcquisition()
savefile = 1;
filename = '20180706_FS_test.txt';
filepath = ['D:\Student_experiments\Virginia\FTSensor\Data\',filename];

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
    targetForce = 5; % [N]
    targetAngles = [0:pi/4:2*pi];
    targetPosx = targetForce*cos(targetAngles)';
    targetPosy = targetForce*sin(targetAngles)';
    rCirCursor = 0.25; % [N]
    rCirTarget = 0.5; % [N]
    
    % Set figure
    hf = figure(1);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %scrsize = get(0, 'Screensize');
    %set(gcf, 'Position', [scrsize(3:4)/4 scrsize(3:4)/2]);
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
    
    % Start acquiring data
    movemtime = 5;
    holdtime = 0.5;
    trialnum = 0;
    cursorholdout = 0;
    countstate = 0;
    state = 'start';
    origstate = 'start';
    global tstart tmovem trelax targetCir tholdstart ht htxt
    h = addlistener(s,'DataAvailable',@(src,event) processData(event,offset,hp));
    s.IsContinuous = true;
    s.Rate = 1000; % scans/sec, samples/sec?
    %s.NotifyWhenDataAvailableExceeds = 100; % Call listener when x samples are available
    s.startBackground();
%     running = true;
%     while running
%         iAngle = randi(8);
%         targetAngle = iAngle;
%         targetCir = circle(rCirTarget,targetPosx(iAngle),targetPosy(iAngle));
%         ht = plot(targetCir(:,1),targetCir(:,2),'r','LineWidth',2);
%         pause(1);
%         delete(ht)
%         quit = input('\Press 0 to stop acquisition.');
%         if ~isempty(quit)
%             running = false;
%         end
%     end
    input('\Press enter to stop acquisition.')
    
    s.stop()
    delete(h)
    close(hf)
else
    disp('No DAQ device detected.')
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
        
%         fc = 200;
%         fs = s.Rate;
%         [b,a] = butter(2,fc/(fs/2),'low');
%         filtData = filter(b,a,event.Data);
        
        offsetData = event.Data-offsetMat;
        
        procData = zeros(size(event.Data));
        for i = 1:size(offsetData,1)
            procData(i,:) = (calibMat*offsetData(i,:)')';
        end
        
        time = event.TimeStamps;
        forceData = procData(:,1:2);
        forceDatax = mean(forceData(:,1));
        forceDatay = mean(forceData(:,2));
        cursorCir = circle(rCirCursor,forceDatax,forceDatay);
        set(hp,'xdata',cursorCir(:,1)','ydata',cursorCir(:,2)');
        
        if strcmp(state,origstate)&&countstate==0
            countstate = countstate+1;
            fprintf('%d ',countstate)
            htxt = text(6.5, 4, state, 'clipping', 'off','Fontsize',12);
        elseif ~strcmp(state,origstate)&&countstate>0
            countstate = 0;
            delete(htxt)
        end
        
        %rectangle('Position',[forceDatax,forceDatay,rCirCursor,rCirCursor],...
        %    'Curvature',[1,1],'FaceColor','g');
        %fill([cursorCir(:,1)]',[cursorCir(:,2)'],'r')

        drawnow;
        
       origstate = state;

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
                    trelax = tic;
                elseif cursorholdout <= 10 && toc(tholdstart)-tmovem > holdtime
                    state = 'success';
                    trelax = tic;
                end
            case 'fail'
                state = 'relax';
                delete(ht)                
            case 'success'
                state = 'relax';
                delete(ht)
            case 'relax'
                if cursorInTarget(cursorCir,circle(rCirCursor+0.1,0,0)) && toc(trelax) > 1
                    state = 'start';
                end
        end

        if savefile
            fprintf(fid, '%d\t %s\t %1.4f\t %1.4f\t %1.4f\n',trialnum,state,time(1),forceDatax,forceDatay);
        end
    end
end
