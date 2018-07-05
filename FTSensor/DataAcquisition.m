% DAQ from F/T Sensor
function DataAcquisition()
savefile = 1;
filename = '20180705_FS_test.txt';
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
    for j = 1:length(targetPosx)
        targetCir = circle(rCirTarget,targetPosx(j),targetPosy(j));
        plot(targetCir(:,1),targetCir(:,2),'r','LineWidth',2)
        %rectangle('Position',[targetPosx(j),targetPosy(j),rCirTarget,rCirTarget],...
        %    'Curvature',[1,1],'FaceColor','r','EdgeColor','r');
    end
    title('2D Force');
    xlabel('F_x [N]'); ylabel('F_y [N]');
    set(gca,'xlim',targetForce*[-1.2 1.2],'ylim',targetForce*[-1.2 1.2]);
    axis square;
    grid on;
    
    % Set file if saving option is enabled
    if savefile
        fid = fopen(filepath,'wt');
        fprintf(fid, ' TimeStamp\t Fx [N]\t Fy [N]\n');
    end
    
    % Start acquiring data
    targetAngle = -1;
    h = addlistener(s,'DataAvailable',@(src,event) processData(event,offset,targetAngle,hp));
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
    function[cursorCir] = processData(event,offset,targetAngle,hp)
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
        set(hp,'xdata',[cursorCir(:,1)]','ydata',[cursorCir(:,2)']);
        %rectangle('Position',[forceDatax,forceDatay,rCirCursor,rCirCursor],...
        %    'Curvature',[1,1],'FaceColor','g');
        %fill([cursorCir(:,1)]',[cursorCir(:,2)'],'r')

        drawnow;
        
        switch state
            case 'start'
                tstart = tic;
                state = 'movement';
            case 'movement'
                if 
        
        if savefile
            fprintf(fid, '%1.4f\t %1.4f\t %1.4f\n',time(1),forceDatax,forceDatay);
        end
    end
end
