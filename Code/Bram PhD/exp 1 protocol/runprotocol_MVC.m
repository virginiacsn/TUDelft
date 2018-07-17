function [] = runprotocol(protocol, settings, startMetTrialNr)

if nargin < 3
    startMetTrialNr = 1;
end

% Init Setup
InitSetup;

% Definieer pointers naar veriabelen in DSPACE
vars = InitAddresses;

% Offsetten
OffSet(vars, settings, protocol);

% Proefpersoon en datum
disp('Protocol instellen.')
% disp(['Krachtniveau: ', num2str(protocol.krachtniveau), ' N'])
protocol.datum = datestr(now, datestr(now,'yyyy-mm-dd HH-MM-SS'));
protocol.proefpersoon = input('Naam van de proefpersoon: ', 's');
protocol.taak='MVC';
  
% Maak data structure
% % Hierin zit zowel het protocol als de gemeten data, zodat deze altijd
% % gekoppeld blijven.
% Ntrials = protocol.settrials;
% totTrials = length(protocol.trials)
% Nsample = settings.meetduur * settings.Fs + 1;
% Force = protocol.krachtniveau;
% Maak lege cell om geheugenruimte voor metingen te reserveren.
% emptycell = cell(Ntrials, 1);
% emptycell(:) = {zeros(Nsample, length(protocol.settrials))};
% Sla protocol en meetdata op in enkele structure.
dataout.protocol = protocol;
% dataout.krachtniveau.volgorde = struct('krachten', zeros(1,length(protocol.krachtniveau)));
% dataout.data = struct('tijd', emptycell, 'FZ', emptycell, 'TX', emptycell,'TY', emptycell);
dataout.data = struct('tijd', emptycell, 'FZ', emptycell, 'TX', emptycell,'TY', emptycell);


voortgang = 0;

% Backup protocol
backupNaam = '.\autobackup\backup_protocol';
save(backupNaam, 'protocol', 'settings');

for trialKracht = 1:length(protocol.krachtniveau)
        mlib('write', vars.boodschap, 'data', 5);
%         dataout.krachtniveau.volgorde(trialKracht) = protocol.krachtniveau(trialKracht);
        pause(2)
for trialNr = startMetTrialNr : Ntrials

    % Laad huidige trial
    currentTrial = protocol.trials(trialNr + ((trialKracht - 1) * protocol.settrials))
    krachtniveau = protocol.krachtniveau(trialKracht);
    
    % Voer trial uit
    [tijd,FZ,TX,TY] = DoeTrial(currentTrial, vars, settings, protocol,krachtniveau);

    % Sla de gemeten data op.
%     dataout.data(trialNr + (21* (trialKracht-1)),trialKracht).tijd(:) = tijd;
%     dataout.data(trialNr + (21* (trialKracht-1)),trialKracht).FZ(:) = FZ;
%     dataout.data(trialNr + (21* (trialKracht-1)),trialKracht).TX(:) = TX;
%     dataout.data(trialNr + (21* (trialKracht-1)),trialKracht).TY(:) = TY;
    
    dataout.data(trialKracht).tijd(:,trialNr) = tijd;
    dataout.data(trialKracht).FZ(:,trialNr) = FZ;
    dataout.data(trialKracht).TX(:,trialNr) = TX;
    dataout.data(trialKracht).TY(:,trialNr) = TY;

    % Sla een reservekopie van dataout op voor eventuele vastlopers of
    % problemen. Met deze data en het protocol is alle meetdata te
    % reconstrueren.
    backupData = dataout.data(trialNr)
    backupNaam = ['.\autobackup\backup_trial_', num2str(trialNr,'%05d')];
    save(backupNaam, 'backupData');
    
    % Bereken voortgang
    voortgang = (trialNr - startMetTrialNr + 1 + (trialKracht - 1) * protocol.settrials) / (totTrials - startMetTrialNr + 1) 
  
%     voortgang = (trialNr - startMetTrialNr + 1) / (Ntrials - startMetTrialNr + 1);
    mlib('write', vars.voortgang, 'data', voortgang)
    
    disp(['Voortgang: ', num2str(trialNr), ' / ', num2str(Ntrials)])

end
end
% Sla de meetdata op
disp('Data wordt opgeslagen')
% saveNaam = ['.\data\', datestr(now,'yyyy-mm-dd HH-MM-SS'), '_', protocol.proefpersoon,'_', num2str(protocol.krachtniveau)];
% saveNaam = ['.\data\', datestr(now,'yyyy-mm-dd HH-MM-SS'), '_', protocol.proefpersoon, '_0_'];
% elseif strcmp(protocol.taak, 'kracht')
% saveNaam = ['.\data\', datestr(now,'yyyy-mm-dd HH-MM-SS'), '_', protocol.proefpersoon, '_999_'];
% else
%     error('taak onduidelijk')
% end
saveNaam = ['.\data\', protocol.proefpersoon, '_', protocol.taak];
saveNaam2 = ['.\data\', protocol.proefpersoon, '_', protocol.taak, '_', 'cont'];
save(saveNaam, 'dataout')
disp('Data is opgeslagen')

% Geef aan dat experiment klaar is
mlib('write', vars.boodschap, 'data', 4)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           INIT SETUP                        %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitSetup
% Load
disp('Loading application to dSpace')
disp(' ')
ok=dos('scoutcmd .\controller\forcereader','-echo');
if ok==0
    disp('Application successfully loaded.')
else
    disp('Error while loading application!')
end
disp(' ')
clear ok

% Init mlib
mlibini
mlib('SelectBoard','ds1005');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           INIT ADDRESSES                      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vars = InitAddresses
vars.FX = mlib('GetTrcVar','Model Root/OffsetsumFX/Out1');
vars.FY = mlib('GetTrcVar','Model Root/OffsetsumFY/Out1');
vars.FZ = mlib('GetTrcVar','Model Root/OffsetsumFZ/Out1');
vars.TX = mlib('GetTrcVar','Model Root/OffsetsumTX/Out1');
vars.TY = mlib('GetTrcVar','Model Root/OffsetsumTY/Out1');
vars.TZ = mlib('GetTrcVar','Model Root/OffsetsumTZ/Out1');
vars.overload = mlib('GetTrcVar','Model Root/Overload warning/Out1');
vars.visueel = mlib('GetTrcVar','Model Root/P:VisSwitch.Gain');
vars.triggerAanUit = mlib('GetTrcVar','Model Root/P:triggerAanUit.Gain');
vars.trigger = mlib('GetTrcVar','Model Root/triggerAanUit/Out1');
vars.resetRecordLight = mlib('GetTrcVar','Model Root/P:reset TriggerSensor.Value');
vars.boodschap = mlib('GetTrcVar','Model Root/P:Boodschap.Value');
vars.voortgang = mlib('GetTrcVar','Model Root/P:Voortgang.Value');
vars.Ftarget = mlib('GetTrcVar', 'Model Root/P:Ftarget.Value');
vars.terug = mlib('getTrcVar','Model Root/P:Terug.Value');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Offset                              %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = OffSet(vars, settings, protocol)

% Maak additionele variabelen aan
add.rawFX = mlib('GetTrcVar','Model Root/Demux/Out1');
add.rawFY = mlib('GetTrcVar','Model Root/Demux/Out2');
add.rawFZ = mlib('GetTrcVar','Model Root/FZgain/Out1');
add.rawTX = mlib('GetTrcVar','Model Root/Demux/Out4');
add.rawTY = mlib('GetTrcVar','Model Root/Demux/Out5');
add.rawTZ = mlib('GetTrcVar','Model Root/Demux/Out6');
add.OffsetFX = mlib('GetTrcVar','Model Root/P:OffsetFX.Value');
add.OffsetFY = mlib('GetTrcVar','Model Root/P:OffsetFY.Value');
add.OffsetFZ = mlib('GetTrcVar','Model Root/P:OffsetFZ.Value');
add.OffsetTX = mlib('GetTrcVar','Model Root/P:OffsetTX.Value');
add.OffsetTY = mlib('GetTrcVar','Model Root/P:OffsetTY.Value');
add.OffsetTZ = mlib('GetTrcVar','Model Root/P:OffsetTZ.Value');



% Zet huidige offsets op nul
mlib('write', add.OffsetFX, 'Data', 0);
mlib('write', add.OffsetFY, 'Data', 0);
mlib('write', add.OffsetFX, 'Data', 0);
mlib('write', add.OffsetTX, 'Data', 0);
mlib('write', add.OffsetTY, 'Data', 0);
mlib('write', add.OffsetTZ, 'Data', 0);


input('Offsetten. Laat de handle los en druk enter.','s');

% Geef de handle tijd om nulpunt te zoeken
disp('Offsetmeting gestart...')
pause(4)

% Duur: 2 seconden
nrDatSamples = 2 * settings.Fs + 1;
downSample = settings.Fds / settings.Fs;

% Initialize MTRACE for data acquisition
mlib('Set',...
    'Trigger','Off',...
    'TraceVars',[add.rawFX; add.rawFY; add.rawFZ; add.rawTX; add.rawTY; add.rawTZ],...
    'AcquisitionMode','single_shot',...
    'DownSampling',downSample,...
    'NumSamples',nrDatSamples);

% start capture (position only) and wait till ready
mlib('StartCapture');

while (mlib('CaptureState')~=0)
    pause(0.01)
end

output = mlib('FetchData');
OffsetFX = mean(output(1,:));
OffsetFY = mean(output(2,:));
OffsetFZ = mean(output(3,:));
OffsetTX = mean(output(4,:));
OffsetTY = mean(output(5,:));
OffsetTZ = mean(output(6,:));

% Stel positie offset waarden in op dspace
mlib('write', add.OffsetFX, 'Data', OffsetFX);
mlib('write', add.OffsetFY, 'Data', OffsetFY);
mlib('write', add.OffsetFZ, 'Data', OffsetFZ);
mlib('write', add.OffsetTX, 'Data', OffsetTX);
mlib('write', add.OffsetTY, 'Data', OffsetTY);
mlib('write', add.OffsetTZ, 'Data', OffsetTZ);

disp('Offsetten: klaar.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Doe enkele trial meting             %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tijd,FZ,TX,TY] = DoeTrial(trial, vars, settings, protocol,krachtniveau)

% Duur van meting instellen
nrDatSamples = settings.meetduur * settings.Fs + 1;
downSample = settings.Fds / settings.Fs;

% % Initialize MTRACE for data acquisition
mlib('Set',...
    'Trigger','Off',...
%     'TriggerVariable', vars.trigger,...
%     'TriggerLevel', 0.1,...
%     'TriggerEdge', 'rising',...
    'TraceVars',[vars.FZ; vars.TX; vars.TY],...
    'AcquisitionMode','continuous',...
    'StreamToDisk','On',...
    'FileName',SaveNaam2,...
    'DownSampling',downSample,...
    'NumSamples',nrDatSamples);
% Initialize MTRACE for data acquisition
% mlib('Set',...
%     'Trigger','On',...
%     'TriggerVariable', vars.trigger,...
%     'TriggerLevel', 0.1,...
%     'TriggerEdge', 'rising',...
%     'TraceVars',[vars.FZ; vars.TX; vars.TY],...
%     'AcquisitionMode','single_shot',...
%     'DownSampling',downSample,...
%     'NumSamples',151);

% Zorgen dat F = 0
disp('Ga naar begin punt (geen kracht)')
mlib('write', vars.boodschap, 'data', 1)
mlib('write', vars.terug, 'data', 1)
Fnow = mlib('read', vars.FZ);
while Fnow > 1
    pause(0.1)
    Fnow = mlib('read', vars.FZ);
end

mlib('write', vars.terug, 'data', 0)
% FIXME: disabled position check
% warning('runprotocol:fixme','Positie-check uitgeschakeld!')
% if strcmp(protocol.taak, 'kracht')
%     pause(1.5)
% else
%     while Xnow > 0
%         pause(0.1)
%         Xnow = mlib('read', vars.X);
%     end
% end

% mlib('write', vars.terug, 'data', 0)

% Stel trialeigenschappen in
% mlib('write', vars.stijfheid, 'data', protocol.stijfheid);
% mlib('write', vars.Fafwijking, 'data', trial.Fafwijking);
mlib('write', vars.visueel, 'data', trial.visueel);
mlib('write', vars.Ftarget, 'data', krachtniveau);


if trial.visueel
    mlib('write', vars.boodschap, 'data', 2);
else
    mlib('write', vars.boodschap, 'data', 3);
end


% Zet trigger op scherp
mlib('write', vars.triggerAanUit, 'data', 1)

% Start getriggerde meting
disp('Meting trial gestart.')
mlib('StartCapture');

while (mlib('CaptureState')~=0)
    pause(0.1)
end

pause(0.5)

% Reset triggersensor (tbv record-lampje)
mlib('write', vars.resetRecordLight, 'data', 1)
pause(0.1)
mlib('write', vars.resetRecordLight, 'data', 0)

% Zet trigger uit en geef aan dat proefpersoon terug moet
mlib('write', vars.triggerAanUit, 'data', 0)
mlib('write', vars.boodschap, 'data', 1)
mlib('write', vars.terug, 'data', 1)

% Lees data uit
output = mlib('FetchData');
disp('Meting trial klaar.')

% X = output(1,:)';
FZ = output(1,:)';
TX = output(2,:)';
TY = output(3,:)';
tijd = (0:nrDatSamples-1)' / settings.Fs;

end
