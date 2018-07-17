function varargout=PosTask_ForcePert_passive(varargin)
global  Achilles Set Get Start Stop Create Remove

% SEE NOTES AT END OF SCRIPT
%% Common: List of Default settings
if nargin==0
    TrialSettings.Description         = 'Relax task - Torque disturbance (MS)';
    %TrialSettings.SignalPeriod        = 43; %keep at this
    TrialSettings.Inertia             = 0.1;
    TrialSettings.Be                  = 0.1;
    TrialSettings.Ke                  = 0;
    TrialSettings.VisTune             = 0.6;
    TrialSettings.ButterCut           = 0.15;
    
    TrialSettings.ScoreThreshold      = 0.5;
    TrialSettings.SignalFreqLo        = 0.2;
    TrialSettings.SignalFreqHi        = 20;
    TrialSettings.SignalFreqTotal     = 20;   % choose if you want to have 40 or 100/200frequencies logarithm. spaced between lo and hi
    TrialSettings.SignalMultiSineGain = 0.8;
    TrialSettings.SignalNTriangles    = 4;
    TrialSettings.SignalTriangleGain  = 0;
    TrialSettings.SignalType          = 1;
    
    TrialSettings.StepTorqueScaling   = 0.1;
    TrialSettings.StepTime            = 3;
    TrialSettings.VisualSigMin        = 0;
    TrialSettings.VisualSigMax        = 0;
    TrialSettings.VisualHorizon       = 3;
    TrialSettings.VisualPast          = 1;
    TrialSettings.VisualType          = 'chirp';
    TrialSettings.VisualSineFreqLo    = 0;
    TrialSettings.VisualSineFreqHi    = 0;
    TrialSettings.VisualBufferSize    = 25;
    TrialSettings.EMG=1;
    varargout={TrialSettings};
end

%% Common: Get protocol & monitor settings out of function input variables
if nargin==3
    TrialSettings=cell2mat(varargin(:,1));
    SessionSettings=cell2mat(varargin(:,2));
    GUISettings=cell2mat(varargin(:,3));
    GUISettings.fcnAddInfo(GUISettings.FeedbackText,[TrialSettings.Description ' Started'])
    
    % Some sessionsetting values are not present
    % ask user if they want to load default values or if they want to
    % import values from a file
    if ~isfield(SessionSettings,'InitialTorque')||~isfield(SessionSettings,'InitialAngle')||~isfield(SessionSettings,'ROM');
        howSvalues = input('Some Session Setting values are not present \n do you want to load default values or import them from a file?\n Type "0" for "default" or "1" for "import":','s');
        if howSvalues=='1'
            [SFileName,SPathName,SFilterIndex] = uigetfile;
            chosenSession = load(fullfile(SPathName,SFileName),'Session');
            SessionSettings  = chosenSession.Session.SessionSettings;
            while (~isfield(SessionSettings,'InitialTorque')||~isfield(SessionSettings,'InitialAngle')||~isfield(SessionSettings,'ROM'))&& howSvalues=='1';
                warning('Chosen an incomplete session, try again')
                howSvalues = input('Some Session Setting values are STILL not present \n do you want to load default values or import them from a file?\n Type "0" for "default" or "1" for "import":','s')
                if howSvalues=='1'
                    [SFileName,SPathName,SFilterIndex] = uigetfile;
                    chosenSession = load(fullfile(SPathName,SFileName),'Session');
                    SessionSettings  = chosenSession.Session.SessionSettings;
                end
            end
            
        else
            % load default values if ROM is not present
            
            if ~isfield(SessionSettings,'ROM');
                ws=Get('workspace');
                curpos=Get('measpos');
                SessionSettings.ROM(1)=ws(1);%max([curpos-0.1 ws(1)]);
                SessionSettings.ROM(2)=ws(2);%min([curpos+0.1 ws(2)]);
                GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Warning: No ROM present, default values loaded')
            end
            
            if ~isfield(SessionSettings,'InitialAngle');
                ws=Get('workspace');
                curpos=Get('measpos');
                SessionSettings.InitialAngle=mean(ws);
                SessionSettings.NeutralAngleEst=mean(ws);
                GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Warning: No zero position present, for Initial Angle and Neutral Angle default values are  loaded')
            end
            
            if ~(isfield(SessionSettings,'MVCP')|isfield(SessionSettings,'MVCD'));
                SessionSettings.MVCP = -25;
                SessionSettings.MVCD = 20;
                GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Warning: No MVCP and/or MVCD present, default values loaded')
            end
            
            if ~isfield(SessionSettings,'InitialTorque');
                SessionSettings.InitialTorque = -1.5;
                GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Warning: No Initial Torque present (gravity), default value is loaded')
            end
            
        end
    end
    %% signal/experiment parameters
    TrialSettings.Ke                  = 50; % aangepast van 20 door Henri/Sven
    TrialSettings.BiasTorque=abs(0.1*SessionSettings.MVCP); 

    VisTune = TrialSettings.VisTune;
    Wn = TrialSettings.ButterCut;
    
    if TrialSettings.SignalMultiSineGain == 0
        Gain = inputdlg('Enter gain of perturbation:','Set gain',1);
        signalgain = str2double(Gain);
        
    else
        signalgain = TrialSettings.SignalMultiSineGain;
    end
    
    if signalgain > 3
        answ = questdlg('Gain > 3: Are you sure?','Warning','Yes','No','No');
        switch answ
            case 'Yes'
                disp('Gain approved')
            case 'No'
                signalgain = 1;
                return
        end
    end
    
    TimeHorizon    = TrialSettings.VisualHorizon;   % seconds
    TimeHist       = TrialSettings.VisualPast;      % seconds history to be seen
    
    % adjust visual trajectory for plantar/dorsal direction
    RangeP=SessionSettings.ROM(1)-SessionSettings.InitialAngle;
    RangeD=SessionSettings.ROM(2)-SessionSettings.InitialAngle;
    
    vmin=min([TrialSettings.VisualSigMin TrialSettings.VisualSigMax]);
    vmax=max([TrialSettings.VisualSigMin TrialSettings.VisualSigMax]);
    na=SessionSettings.InitialAngle;
    
    if vmin>=0
        Ftrajmin = na+vmin*RangeD;
        Ftrajmax = na+vmax*RangeD;
    elseif vmax<=0
        Ftrajmin = na-vmin*RangeP;
        Ftrajmax = na-vmax*RangeP;
    elseif vmin*vmax<0;
        Ftrajmin = na-vmin*RangeP;
        Ftrajmax = na+vmax*RangeD;
    end
    
    buffersize= TrialSettings.VisualBufferSize ;
    
    if strcmp(TrialSettings.VisualType,'chirp')
        FreqTraj =  [TrialSettings.VisualSineFreqLo TrialSettings.VisualSineFreqHi] ;
    else
        FreqTraj =  [TrialSettings.VisualSineFreqLo];
    end
    
    Ie = TrialSettings.Inertia;   %environment inertia
    Be = TrialSettings.Be; %environment damping coefficient
    Ke = TrialSettings.Ke; %environment stiffness
    
    % Load perturbation signal
    filepath = 'D:\Student_experiments\Achilles - BEP (Mark)\AchFuncs\Bo\40_frequencies\';
    filey1 = 'finished_Hz40_40points.mat';
    Filey1 = fullfile(filepath,filey1);
    load(Filey1,'Multisine');
    y0 =zeros(1,2048*5);
    y1 =  Multisine'.*signalgain;
    y1 = [y0 y1 y0];
    
    save('torque','y1');
    Drivefile = 'torque.mat';%['Signals\',SigNames{str2double(signr)}];
    Drivefilegain = 1;
    Drivefilefadingduration = 0;
    T=(length(y1))/2048;
    
    if GUISettings.offline;T=10;end
    
    %% Initialise Porti7 & Achilles
    
    EMGtry                                                      % loads TMSi code and initializes 'Device'
    GUISettings.fcnAddInfo(GUISettings.FeedbackText,messageEMG) % feedback from EMG-initialization
    Startposition=SessionSettings.InitialAngle;
    
    Remove('all')
    Create('spring','MySpring');
    Create('damper','MyDamper');
    Create('biastorque','MyBias')
    Set('torquecompensation','disable')
    
     GUISettings.fcnSetPosCos(Startposition,2,GUISettings.fcnAddInfo,GUISettings.FeedbackText);

    Set('MySpring','pos',Startposition,'stiffness',Ke,'effectmaxtorque',100,'dampfactor',0);
    Set('MyDamper','dampcoef',Be,'coulombthreshold',10);
    Set('inertia',Ie);
    Set('maxtorque',150);
    Set('MyBias','torque',TrialSettings.BiasTorque)
    Start('MySpring','MyDamper','MyBias');
    
    % Drivefile settings
    Set('drivefile',Drivefile,'torque',...
        'drivefile gain',Drivefilegain,...
        'drivefile fadingduration',Drivefilefadingduration,...
        'state','torque');
    
    % DataloggerFlush Settings
    Set('dataloggerFlush',512,1,...
        'sampletime',...
        '#Achilles.DriveFile.Output torque',...
        '#Achilles.Analogue IO.Analogue_1.Constant value',...
        '#Achilles.VelCmdScaling.input',...
        '#Achilles.Meas position',...
        '#Achilles.Meas velocity',...
        '#Achilles.Meas torque',...
        '#Achilles.Meas torque compensated',...
        '#Achilles.DriveFile.ReductionFactorFader.FaderValue');
    
    % dataloggerThread settings
    Set('dataloggerThread',2,TrialSettings.logname,...
        'sampletime',...
        '#Achilles.DriveFile.Output torque',...
        '#Achilles.Analogue IO.Analogue_1.Constant value',...
        '#Achilles.VelCmdScaling.input',...
        '#Achilles.Meas position',...
        '#Achilles.Meas velocity',...
        '#Achilles.Meas torque',...
        '#Achilles.Meas torque compensated',...
        '#Achilles.DriveFile.ReductionFactorFader.FaderValue');
    
    %% Visualization initialization
    
    t_prep=1;
    t_hold=1;
    fsv=500; %sample rate for visualization
    dotsize=0.08;
    ndots=20;
    pheadsize=0.25;
    peyesize=0.05;
    Cdots=[1;0]*linspace(-1,1,ndots);
    t_traj=(0:1/fsv:T+t_prep+t_hold).';
    CursorPos=-1+2*TimeHist/TimeHorizon;
    hcf=gcf;
    pos_vis=0;
    [subjecttext,sth,hphead,hpeyes,hpmouth,hdots]=InitFigure(GUISettings.ParentAxesSubject,GUISettings.AxesSubject,CursorPos,Cdots,pheadsize,peyesize,dotsize,ndots);
    set(subjecttext,'String','Volg de bolletjes...')
    
    Nsteps=1;Nreps=1;UpDown=1;
    
    T_signal=MakeTrajTV(TrialSettings.VisualType,Ftrajmin,Ftrajmax,FreqTraj,Nsteps,Nreps,UpDown,t_traj(1:T*fsv));
    T_traj=[T_signal(1)*ones(t_prep*fsv,1);T_signal;T_signal(end)*ones(t_hold*fsv,1)];
    %pvis=polyfit([SessionSettings.AROM(1),SessionSettings.AROM(2)],[-0.8 0.8],1);
    dif1 = na -SessionSettings.ROM(1);
    dif2 = SessionSettings.ROM(2)-na;
    range1 = na-VisTune*dif1;
    range2 = na+VisTune*dif2;
    pvis=polyfit([range1,range2],[-0.8 0.8],1);
    %drawnow;
    
    %% Start logger
    
    pmsound1 = audioread([cd '\Sounds\pacman_intro2.wav']);
    ap=audioplayer(pmsound1,11025);
    play(ap);
    
    % Monitor
    if GUISettings.MonitorOn==1;
        set(GUISettings.ButToggleMonitor,'Value',1);
    end
    
    set(GUISettings.ButExit,'String','Stop')
    
    if get(GUISettings.ButToggleMonitor,'Value')==1
        GUISettings.looppause=0;
        Set('dataloggerFlush',512,4,'sampletime','#Achilles.Meas position','#Achilles.Meas torque','#Achilles.DriveFile.Output torque');
        Achilles.matrixColumnCount=4; %% needed to prevent growth of data matrix->contact MOOG, something with workspaces and functions
        Start('dataloggerFlush');
    end
    
    if EMGEnabled
        sampler.start()
    else
        EMGdatout=[];
    end
    
    %% recording loop
    
    loggerstarted=0;
    pos_buffer=zeros(1,buffersize);
    v_buffer=zeros(1,buffersize);
    k_buffer=zeros(1,buffersize);
    [bfilp,afilp]=butter(2,Wn);
    [bfilv,afilv]=butter(2,Wn);
    Ttot=T+t_prep+t_hold;
    score=0;
    
    Start('dataloggerThread','drivefile');
    GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Recording...')
    ws=Get('workspace');
    tic
    
    tt = 0; tt2=0;
    try
        HA_open(Achilles.AchillesIPAddress)
        %HA_command('set digital_output_1 true');
        HA_command('set analogue_io.analogue_1.enabled false');
        %HA_command('get analogue_io.analogue_1.enabled');
        HA_command('set analogue_io.analogue_1.output_signal_type constant');
        %HA_command('get analogue_io.analogue_1.output_signal_type');
        HA_command('set analogue_io.analogue_1.constant_value -1');
        %HA_command('get analogue_io.analogue_1.constant_value');
        
        while toc<Ttot+1 && get(GUISettings.ButStop,'value')==0;
            %Set dig out1 to HIGH -> trigger
            if (tt == 0) && (toc > (Ttot-38))
                HA_command('set analogue_io.analogue_1.enabled true');
                HA_command('set analogue_io.analogue_1.constant_value 5');
                %HA_command('set digital_output_1 false');
                tt = tt + 1;
            end
            
            % update waitbar in protocol player
            dum=get(GUISettings.Progressbar,'Position');
            dum(3)=toc/Ttot;
            set(GUISettings.Progressbar,'Position',dum);
            % update signals if signal monitor is switched on
            if get(GUISettings.ButToggleMonitor,'Value')==1
                GUISettings.PlotBuffer=fcnRefreshMonitor(GUISettings);
            end
            % update EMG signals
            if EMGEnabled
                samples = sampler.sample();
                emg_data.append(samples(channel_subset,:));
            end
            
            % update visualization
            pos_meas=Get('measpos');
            %             pos_meas=getinput(hcf,pos_filt(end));
            pos_prev=pos_vis;
            pos_buffer=[pos_buffer(2:end) pos_meas];
            %             pos_filt=filter(bfilp,afilp,pos_buffer);
            pos_filt=filter(bfilp,afilp,pos_buffer-pos_buffer(1))+pos_buffer(1);
            
            pos_now=pos_filt(end);
            
            pos_vis=polyval(pvis,pos_now);%pos_now;%2*(pos_now-ws(1))/(ws(2)-ws(1))-1;
            
            v_meas=pos_vis-pos_prev;
            v_buffer=[v_buffer(2:end) v_meas];
            %             v_filt=filter(bfilv,afilv,v_buffer);
            v_filt=filter(bfilp,afilp,v_buffer-v_buffer(1))+v_buffer(1);
            
            v_now=v_filt(end);
            
            %             figure(8)
            %             subplot(211)
            %             plot(1:buffersize,pos_filt,1:buffersize,pos_buffer)
            %             ylim([-1 1])
            %             subplot(212)
            %             plot(1:buffersize,v_filt,1:buffersize,v_buffer)
            %             k_elaps=round(toc*fsv)+1;
            k_elaps=min([(T+t_prep+t_hold)*fsv round(toc*fsv)+1]);
            k_buffer=[k_buffer(2:end) k_elaps];
            
            theta=atan(v_now/(range(k_buffer)/(TimeHorizon*fsv)));
            R=[cos(theta) -sin(theta);sin(theta) cos(theta)];
            if k_elaps>t_prep*fsv&&loggerstarted==0
                set(subjecttext,'Visible','off')
                
                loggerstarted=1;
                
            end
            
            %% Refresh visualization
            TrajBusyStartIdx=max([1 k_elaps-TimeHist*fsv+1]);
            TrajBusyEndIdx=min([length(T_traj) k_elaps+(TimeHorizon-TimeHist)*fsv]);
            
            T_trajBusy=T_traj(TrajBusyStartIdx:TrajBusyEndIdx);
            if TrajBusyEndIdx==length(T_traj)
                T_trajEnd=T_traj(end)*ones(TimeHorizon*fsv-size(T_trajBusy,1),1);
            else
                T_trajEnd=[];
            end
            if TrajBusyStartIdx==1
                T_trajStart=T_traj(1)*ones(TimeHorizon*fsv-size(T_trajBusy,1),1);
            else
                T_trajStart=[];
            end
            
            ya=([T_trajStart;T_trajBusy;T_trajEnd])';
            yd=polyval(pvis,ya);
            
            Cdots(1,:)=1-2*rem(rem(k_elaps,(TimeHorizon)*fsv)/((TimeHorizon)*fsv)+linspace(0,1,ndots),1);
            dum=(Cdots(1,:)+1)/2*(length(yd)-1);
            Cdots(2,:)=(dum-floor(dum)).*yd(ceil(dum)+1)+(ceil(dum)-dum).*yd(floor(dum)+1);
            
            sensx=.25;
            sensy=.4;
            for ii=1:ndots
                set(hdots(ii),'Position',[Cdots(1,ii)-.5*dotsize Cdots(2,ii)-.5*dotsize dotsize dotsize])
                if Cdots(1,ii)>CursorPos-sensx*pheadsize && Cdots(1,ii)<CursorPos+sensx*pheadsize && Cdots(2,ii)<pos_vis+sensy*pheadsize && Cdots(2,ii)>pos_vis-sensy*pheadsize
                    if strcmp(get(hdots(ii),'Visible'),'on')
                        score=str2double(get(sth,'String'));
                        if toc>t_prep && toc <(T+t_prep+t_hold)
                            set(sth,'String',num2str(score+1));
                        end
                    end
                    set(hdots(ii),'Visible','off');
                    %                     play(ap2);
                end
                if Cdots(1,ii)>0.95
                    set(hdots(ii),'Visible','on');
                end
            end
            
            set(hphead,'Position',[CursorPos-.5*pheadsize pos_vis-.5*pheadsize pheadsize pheadsize])
            xydum=R*[[.5*pheadsize .5*pheadsize];[.5*pheadsize -.5*pheadsize]];
            eyedum=R*[-.5*peyesize;.2*pheadsize];
            set(hpeyes,'Position',[ [CursorPos pos_vis]+eyedum' peyesize peyesize]);
            set(hpmouth,'XData',[CursorPos CursorPos+xydum(1,:)],'YData',[pos_vis pos_vis+xydum(2,:)])
            drawnow;
            
            %Set dig out1 to LOW -> trigger
            if (tt2 == 0) && (toc > (Ttot-8))
                HA_command('set analogue_io.analogue_1.constant_value -1');
                HA_command('set analogue_io.analogue_1.enabled false');
                %HA_command('set digital_output_1 true');
                tt2 = tt2 + 1;
            end
        end
        HA_close()
    catch matlabException
        warning(sprintf('Error while trying to execute the recording loop. \n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
        message = '';
        HA_command('set analogue_io.analogue_1.constant_value -1');
        HA_command('set analogue_io.analogue_1.enabled false');
        %HA_command('set digital_output_1 true');
        HA_close()
    end
    if get(GUISettings.ButStop,'value')==1;
        set(GUISettings.ButStop,'value',0);
        GUISettings.fcnAddInfo(GUISettings.FeedbackText,'Stop button pressed, recording interupted')
    end
    
    %% Finish and store data
    set(GUISettings.ButToggleMonitor,'Value',0);
    if GUISettings.offline
        TrialSettings.logname='c:\temp\dumlog';
        Stop('dataloggerThread','drivefile');
        load(TrialSettings.logname);
        M=randn(1,10);
    else
        Stop('dataloggerThread','drivefile');
        M=load(TrialSettings.logname);
        
        OutputInfo.DataRecorded=1;
    end
    
    if EMGEnabled
        % Step 11: Stop sampler.
        sampler.stop();
        % Step 12: Disconnect with device.
        sampler.disconnect();
        emg_data.trim();
    end
    
    %% EMG for plotting
    if EMGEnabled
        emgpp = emg_data.samples;
        %         emgpp=fcnPostProcEMG(Channels,EMGdatout,DeviceSampleRate);
    else
        emgpp=zeros(4,1);
    end
    %% outputs:
    
    OutputInfo.EMGEnabled=EMGEnabled;
    OutputInfo.offline=GUISettings.offline;
    OutputInfo.emgpp=emgpp;
    OutputInfo.EMGtotal = emg_data;
    score;
    TotDots=ndots*T/TimeHorizon+ndots*TimeHist/TimeHorizon;
    
    OutputData.AchData=M;
    OutputData.EMGData=emg_data.samples;
    OutputData.Drivefile=y1;
    OutputData.TrialSettings=TrialSettings;
    OutputData.Score=score;
%     Standardeviation = std(OutputData.AchData(1024*14:1024*43,9).*(180/pi))
    varargout={OutputInfo,OutputData,SessionSettings};
    
    % Step 13: Cleanup library.
    library.destroy();
    
end
end

function [subjecttext,sth,hphead,hpeyes,hpmouth,hdots]=InitFigure(ParentAxesSubject,AxesSubject,CursorPos,Cdots,pheadsize,peyesize,dotsize,ndots)

%     pheadsize=0.5;
%     peyesize=0.1;

if isfield(get(ParentAxesSubject),'Name')
    figure(ParentAxesSubject);
    cla
    haxes=AxesSubject;
else
    haxes=AxesSubject;
    cla(haxes);
end
set(haxes,'Color',[0 0 0])

%%%Position of the figure
%Get screen positions
set(0,'Units','pixels'); scnpos = get(0,'MonitorPosition'); %get screen positions (NOTE: situation at startup of Matlab!!!, so if changed later -> restart Matlab!)
if scnpos(1,1) ~= 1
    scnpos = [scnpos(2,:); scnpos(1,:)]; %swap M1 and M2, because M1 is not the first sometimes??? -> don't know why...
end
mon_nr = size(scnpos, 1); %number of monitors

%Monitor 1 (Main screen)
y_offset = 0;       %vertical offset [pixels] of positionvector (bottom of figure)
y_taskbar = 40;     %vertical size taskbar [pixels] (only used for screen 1)

x1 = scnpos(1,1);
y1 = scnpos(1,2) + y_taskbar + y_offset;
fsize_x = scnpos(1,3);
fsize_y = (scnpos(1,4) - y_taskbar - y_offset);
scnpos1_real = [x1,     y1,     fsize_x,   fsize_y]; %define usable the size and location of the figures

if mon_nr > 1
    %Monitor 2 (Second screen -> Visualization experiment)
    x1 = scnpos(2,1);
    y1 = scnpos(2,2);
    fsize_x = scnpos(2,3);
    fsize_y = scnpos(2,4);
    scnpos2_real = [x1,     y1,     fsize_x,   fsize_y]; %define the size and location of the figures
else
    scnpos2_real = scnpos1_real;
    disp('Only one monitor connected!!')
end
%set(hfig,'Outerposition',scnpos1_real)
set(ParentAxesSubject,'units','pixels','Position',scnpos2_real)

hphead = rectangle('Parent',haxes,...
    'Position',[CursorPos-.5*pheadsize -.5*pheadsize pheadsize pheadsize],...
    'Curvature',[1 1],...
    'FaceColor','y');

hpmouth=fill([CursorPos CursorPos+.5*sqrt(2)*pheadsize CursorPos+.5*sqrt(2)*pheadsize],...
    [0 .5*sqrt(2)*pheadsize -.5*sqrt(2)*pheadsize],...
    [0 0 0],...
    'Parent',haxes,...
    'linewidth',1,...
    'edgecolor',[0 0 0]);
if isfield(get(ParentAxesSubject),'Name')
    set(ParentAxesSubject,'Name','Achilles View',...
        'MenuBar','none',...
        'NumberTitle','off')
end
hdots=zeros(ndots,1);
for ii=1:ndots;
    hdots(ii)=rectangle('Parent',haxes,...
        'Position',[Cdots(1,ii)-.5*dotsize Cdots(2,ii)-.5*dotsize dotsize dotsize],...
        'Curvature',[1 1],...
        'FaceColor','r');
end



hpeyes = rectangle('Parent',haxes,...
    'Position',[CursorPos-.5*peyesize .2*pheadsize peyesize peyesize],...
    'Curvature',[1 1],...
    'FaceColor','k');

sth = uicontrol('Parent',ParentAxesSubject,...
    'Style','text',...
    'Units','Normalized',...
    'String','0',...
    'Position',[0 0.85 0.2 0.1],...
    'FontSize',20,...
    'BackGroundColor',[0 0 0],...
    'ForegroundColor',[1 1 1],...
    'Visible','on');
subjecttext = uicontrol('Parent',ParentAxesSubject,...
    'Style','text',...
    'Units','Normalized',...
    'String','...',...
    'Position',[0.25 0.1 0.5 0.2],...
    'Visible','on',...
    'BackgroundColor',[0 0 0],...
    'ForegroundColor',[1 1 1],...
    'FontSize',24);

end

function pos_input=getinput(hcf,pos_input)
num=double(get(hcf,'CurrentCharacter'));
if num==31
    pos_input=pos_input-0.1;
elseif num==30
    pos_input=pos_input+0.1;
end
set(hcf,'CurrentCharacter','o')
end


%% NOTES
% The analogue output is used as a trigger instead of the digital
% output. Reason: digital output when connected to Porti 7 gives a
% 1.5V output instead of 5V -> trigger is not recorded on Porti 7.
% It takes time before the value of the trigger output changes in
% values. -> it takes time before the EMG is triggerd
% -> make sure that the first and last second of the experimental
% protocol do not contain a crucial part for which you want to
% record the EMG
%
% It is important to always use HA_close after HA_open.