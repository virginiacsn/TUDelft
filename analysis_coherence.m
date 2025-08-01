%% COHERENCE ANALYSIS
Aparams.cohparams.data = 'app';
Aparams.cohparams.tseg = 1;
Aparams.cohparams.nseg = 10;
Aparams.cohparams.my_nseg = 10;
Aparams.cohparams.window = @(N) hanning(N);
Aparams.cohparams.CLoverlap = 1;
fields_coh = {'rect'};
fc = 80;

if strcmp(Aparams.cohparams.data,'avg')
    trial_coh_force = cohStruct(trial_avg_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_coh_EMG = cohStruct(trial_avg_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
else
    trial_coh_force = cohStruct(trial_app_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_coh_EMG = cohStruct(trial_app_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
end

%% Field and channels to plot
% Number of muscles
nmusc = length(EMGparams.channelName)-1;
% Number of muscle pair combinations
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2

field = 'rect';
freqBand = {'Alpha','Beta','Gamma'};

% Find muscle in all muscle pair combinations
plotBB = sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],'BB'),[2,nmusccomb]));
plotTLH = sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],'TLH'),[2,nmusccomb]));
plotDA = sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],'DA'),[2,nmusccomb]));
plotDP = sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],'DP'),[2,nmusccomb]));

% Find indexes of muscle pair combinations to plot. Combinations of control
% muscles.
plotmusc = find((plotBB&plotDA)|(plotDP&plotBB)|(plotDA&plotTLH));%|(plotDA&plotDP)|(plotBB&plotTLH)|(plotTLH&plotDP));

% Coherence (my_coh) limits for each field in stuct. coh_lim(i,j) will be max coh 
% level between force-control and EMG-control where rows are target angles 
% to compare and cols muscle combinations (all)
for h = 1:length(fields_coh)
    for i = 1:length(Aparams.angComp)
        for j = 1:nmusccomb
            coh_lim.(fields_coh{h})(i,j) = max([max(trial_coh_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).(fields_coh{h}).my_coh(3:end-3,j)),...
                max(trial_coh_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).(fields_coh{h}).my_coh(3:end-3,j))]);
        end
    end
end

%% Coherence 
%% Fig per muscle pair in plotmusc. Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and targets from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence: Musc: ',trial_coh_force(1).(field).muscles{plotmusc(j)}{1},...
        ',',trial_coh_force(1).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_coh_force(iangf).(field).my_fcoh(:,plotmusc(j)),trial_coh_force(iangf).(field).my_coh(:,plotmusc(j)));
        hold on;
        line(xlim,trial_coh_force(iangf).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_coh_EMG(iangE).(field).my_fcoh(:,plotmusc(j)),trial_coh_EMG(iangE).(field).my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_coh_EMG(iangE).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(i,plotmusc(j))]);
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair in plotmuscpair. Subplot per target from angCompPair. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and
% targets from angCompPair
for i = 1:length(Aparams.angCompPair)
    plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(1)),...
        [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(3)),[2,nmusccomb])));
    
    figure('Name',['My Coherence: Musc: ',trial_coh_force(1).(field).muscles{plotmuscpair}{1},...
        ',',trial_coh_force(1).(field).muscles{plotmuscpair}{2}]);
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    for k = 1:length(Aparams.angCompPair{i})
        iangf = find([trial_coh_force.angle] == Aparams.angCompPair{i}(k));
        iangE = find([trial_coh_EMG.angle] == Aparams.angCompPair{i}(k));
        
        if rem(length(Aparams.angCompPair{i}),2) == 0
            subplot(length(Aparams.angCompPair{i})/2,2,k);
        else
            subplot(length(Aparams.angCompPair{i}),1,k);
        end
        h1 = plot(trial_coh_force(iangf).(field).my_fcoh(:,plotmuscpair),trial_coh_force(iangf).(field).my_coh(:,plotmuscpair),'b','Linewidth',2);
        hold on;
        line(xlim,trial_coh_force(iangf).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_coh_EMG(iangE).(field).my_fcoh(:,plotmuscpair),trial_coh_EMG(iangE).(field).my_coh(:,plotmuscpair),'r','Linewidth',2);
        hold on;
        line(xlim,trial_coh_EMG(iangE).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(Aparams.angCompUni==Aparams.angCompPair{i}(k),plotmuscpair)]);
        if k == length(Aparams.angCompPair{i})
            xlabel('Frequency [Hz]');
        end
        ylabel('Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angCompPair{i}(k))),' (',Aparams.muscCompPair{i}{k},')']);
        if k == 1
            legend([h1,h2],'FC','EC')
        end
        set(gca,'FontSize',12);
    end
end

%% Before in one fig
figure('Name','Coherence');
h = 0;
sangpair = {Aparams.angCompPair{2},Aparams.angCompPair{3},Aparams.angCompPair{1}};
smuscpair = {Aparams.muscCompPair{2},Aparams.muscCompPair{3},Aparams.muscCompPair{1}};

for k = 1:length(Aparams.angCompPair)
    for i = 1:length(Aparams.angCompPair{i})
        
        plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(1)),...
            [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(3)),[2,nmusccomb])));
        
        iangf = find([trial_coh_force.angle] == sangpair{i}(k));
        iangE = find([trial_coh_EMG.angle] == sangpair{i}(k));
        h = h+1;
        subplot(length(Aparams.angCompPair),length(Aparams.angCompPair),h);
        
        h1 = plot(trial_coh_force(iangf).(field).my_fcoh(:,plotmuscpair),trial_coh_force(iangf).(field).my_coh(:,plotmuscpair),'b','Linewidth',2);
        hold on;
        line(xlim,trial_coh_force(iangf).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_coh_EMG(iangE).(field).my_fcoh(:,plotmuscpair),trial_coh_EMG(iangE).(field).my_coh(:,plotmuscpair),'r','Linewidth',2);
        hold on;
        line(xlim,trial_coh_EMG(iangE).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(Aparams.angCompUni==sangpair{i}(k),plotmuscpair)]);
        if k == length(Aparams.angCompPair{i})
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel('Coh [-]');
        end
        title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']);
        if i == length(Aparams.angCompPair{i}) && k == 1
            legend([h1,h2],'FC','EC')
        end
        set(gca,'FontSize',13);
    end
end

%% Fig per muscle pair (all). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (all) between tasks and targets from angComp
for j = 1:nmusccomb
    figure('Name','My Coherence');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_coh_force(iangf).(field).my_fcoh(:,j),trial_coh_force(iangf).(field).my_coh(:,j));
        hold on;
        line(xlim,trial_coh_force(iangf).(field).my_CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_coh_EMG(iangE).(field).my_fcoh(:,j),trial_coh_EMG(iangE).(field).my_coh(:,j),'r');
        hold on;
        line(xlim,trial_coh_EMG(iangE).(field).my_CL(j)*[1 1],'Color','k','LineStyle','-.');
                
        xlim([trial_coh_force(iangf).(field).my_fcoh(3,j) fc]);
        %ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (plotmusc). Subplot per task (col) per target from angComp (row). Line per processing method.
% Compare processing method (filt/rect) for each muscle combination between
% task and target from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence - Comparison Proc; Musc: ',trial_coh_force(1).filt.muscles{plotmusc(j)}{1},...
        ',',trial_coh_force(1).filt.muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    h = 1;
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        subplot(length(Aparams.angComp),2,h);
        h1 = plot(trial_coh_force(iangf).filt.my_fcoh(:,plotmusc(j)),trial_coh_force(iangf).filt.my_coh(:,plotmusc(j)),'b');
        hold on;
        line(xlim,trial_coh_force(iangf).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_coh_force(iangf).rect.my_fcoh(:,plotmusc(j)),trial_coh_force(iangf).rect.my_coh(:,plotmusc(j)),'c');
        line(xlim,trial_coh_force(iangf).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max([coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j))])]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        subplot(length(Aparams.angComp),2,h+1);
        h1 = plot(trial_coh_EMG(iangE).filt.my_fcoh(:,plotmusc(j)),trial_coh_EMG(iangE).filt.my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_coh_EMG(iangE).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_coh_EMG(iangE).rect.my_fcoh(:,plotmusc(j)),trial_coh_EMG(iangE).rect.my_coh(:,plotmusc(j)),'m');
        line(xlim,trial_coh_EMG(iangE).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max([coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j))])]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        h = h+2;
    end
end

%% Significant coherence
%% Fig of area of coherence per muscle pair (plotmuscpair). Polar plot.
for i = 1:length(Aparams.angCompPair)
    plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(1)),...
        [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(3)),[2,nmusccomb])));
    
    figure('Name',['Area of significant coherence: Musc: ',trial_coh_force(1).(field).muscles{plotmuscpair}{1},...
        ',',trial_coh_force(1).(field).muscles{plotmuscpair}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    fsig = zeros(length(Aparams.angCompPair{i}),3);
    Esig = zeros(length(Aparams.angCompPair{i}),3);
    for k = 1:length(Aparams.angCompPair{i})
        iangf = find([trial_coh_force.angle] == Aparams.angCompPair{i}(k));
        iangE = find([trial_coh_EMG.angle] == Aparams.angCompPair{i}(k));
        
        fsig(k,:) = trial_coh_force(iangf).rect.asig_coh(1:3,plotmuscpair)';
        Esig(k,:) = trial_coh_EMG(iangE).rect.asig_coh(1:3,plotmuscpair)';
    end
    
    for b = 1:length(freqBand)
        subplot(1,length(freqBand),b)
        polarscatter(Aparams.angCompPair{i},fsig(:,b),60,'filled')
        hold on
        polarscatter(Aparams.angCompPair{i},Esig(:,b),60,'filled','r')
        thetaticks(rad2deg(Aparams.angCompPair{i})); thetaticklabels(rad2deg(Aparams.angCompPair{i}));
        set(gca,'FontSize',12);
        if b == 3
            legend('ForceCO','EMGCO','Location','bestoutside');
        end
        title(freqBand{b});
    end
end

%% Fig of area of coherence for muscle pair (plotmuscpair). Polar plot.
figure('Name','Area of significant coherence');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h = 0;
for i = 1:length(Aparams.muscCompPair)
    plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(1)),...
        [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(3)),[2,nmusccomb])));
    
    fsig = [];
    Esig = [];
    
    for k = 1:length(Aparams.angCompPair{i})
        iangf = find([trial_coh_force.angle] == Aparams.angCompPair{i}(k));
        iangE = find([trial_coh_EMG.angle] == Aparams.angCompPair{i}(k));
        
        fsig(k,:) = trial_coh_force(iangf).(field).asig_coh(1:3,plotmuscpair)';
        Esig(k,:) = trial_coh_EMG(iangE).(field).asig_coh(1:3,plotmuscpair)';
    end
    
    for j = 1:length(freqBand)
        h = h+1;
        subplot(length(Aparams.angCompPair),length(freqBand),h)
        g = polar(0,max(max(abs([fsig;Esig]))),'w');
        polarticks(8,[g]);       
        hold on;
        f = polar((Aparams.angCompPair{i}),fsig(:,k)','b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar((Aparams.angCompPair{i}),Esig(:,k)','r-.o');
        e.MarkerFaceColor = 'r';      
        if i == 1
            title(freqBand{j})
        end
        if j == 1
            ylabel([Aparams.muscCompPair{i}{1},'-',Aparams.muscCompPair{i}{3}],'FontWeight','bold');
        end
        set(gca,'FontSize',14);
    end
end

%% Fig of area of coherence for muscle pair (plotmuscpair). Plot.
figure('Name','Area of significant coherence');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h = 0;
for i = 1:length(Aparams.muscCompPair)
    plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(1)),...
        [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],Aparams.muscCompPair{i}(3)),[2,nmusccomb])));
    
    fsig = [];
    Esig = [];
    
    for k = 1:length(Aparams.angCompPair{i})
        iangf = find([trial_coh_force.angle] == Aparams.angCompPair{i}(k));
        iangE = find([trial_coh_EMG.angle] == Aparams.angCompPair{i}(k));
        
        fsig(k,:) = trial_coh_force(iangf).(field).asig_coh(1:3,plotmuscpair)';
        Esig(k,:) = trial_coh_EMG(iangE).(field).asig_coh(1:3,plotmuscpair)';
    end
    
    for j = 1:length(freqBand)
        h = h+1;
        subplot(length(Aparams.angCompPair),length(freqBand),h)
        f = plot(rad2deg(Aparams.angCompPair{i}),fsig(:,j),'b-.o','MarkerFaceColor','b');
        hold on;
        e = plot(rad2deg(Aparams.angCompPair{i}),Esig(:,j),'r-.o','MarkerFaceColor','r');
        if i == 1
            title(freqBand{j})
        elseif i == length(Aparams.muscCompPair)
            xlabel('Target [deg]');
        end
        if j == 1
            ylabel([Aparams.muscCompPair{i}{1},'-',Aparams.muscCompPair{i}{3}],'FontWeight','bold');
        end
        if j == 3 && i == 1
            legend([f,e],'FC','MC');
        end
        set(gca,'XTick',rad2deg(Aparams.angCompPair{i}));
        set(gca,'FontSize',13);
    end
end

%% Fig of normalized area of coherence for muscle pair (plotmuscpair). Plot.
sangpair = {Aparams.angCompPair{2},Aparams.angCompPair{3},Aparams.angCompPair{1}};
smuscpair = {Aparams.muscCompPair{2},Aparams.muscCompPair{3},Aparams.muscCompPair{1}};

figure('Name','Normalized area of significant coherence');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h = 0;
for j = 1:length(freqBand)
    for i = 1:length(Aparams.muscCompPair)
        plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(1)),...
            [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(3)),[2,nmusccomb])));        
        fsig = [];
        Esig = [];
        
        for k = 1:length(Aparams.angCompPair{i})
            iangf = find([trial_coh_force.angle] == sangpair{i}(k));
            iangE = find([trial_coh_EMG.angle] == sangpair{i}(k));
            
            fsig(k,:) = trial_coh_force(iangf).(field).nasig_coh(1:3,plotmuscpair)';
            Esig(k,:) = trial_coh_EMG(iangE).(field).nasig_coh(1:3,plotmuscpair)';
        end        
        h = h+1;
        subplot(length(Aparams.angCompPair),length(freqBand),h)
        f = plot(rad2deg(sangpair{i}),fsig(:,j),'b-.o','MarkerFaceColor','b');
        hold on;
        e = plot(rad2deg(sangpair{i}),Esig(:,j),'r-.o','MarkerFaceColor','r');
        if j == 1
            title([smuscpair{i}{1},'-',smuscpair{i}{3}],'FontWeight','bold');          
        elseif j == length(Aparams.muscCompPair)
            xlabel('Target [deg]');
        end
        if i == 1
            ylabel({['\bf',freqBand{j}];'\rmArea Coh [-]'},'interpreter','tex');
        end
        if j == 1 && i == 3
            legend([f,e],'FC','MC');
        end
        xlim(rad2deg([sangpair{i}(1)-0.2 sangpair{i}(3)+0.2]));
        set(gca,'XTick',rad2deg(sangpair{i}));
        set(gca,'FontSize',13);
    end
end

