function [] = runprotocol(protocol, settings, startMetTrialNr)

if nargin < 3
    startMetTrialNr = 1;
end

% Init Proprio
InitProprio;

% Definieer pointers naar veriabelen in DSPACE
vars = InitAddresses;

% Offsetten
OffSet(vars, settings);

% Proefpersoon en datum
disp('Protocol instellen.')
disp(['Stijfheid: ', num2str(protocol.stijfheid), ' N/m'])
protocol.datum = datestr(now, datestr(now,'yyyy-mm-dd HH-MM-SS'));
protocol.proefpersoon = input('Naam van de proefpersoon: ', 's');
taakinvoer = input('Taakinstructie: krachttaak/positietaak [k/p]', 's');
if taakinvoer=='k'
    mlib('write', vars.isKrachttaak, 'Data', 1);
    protocol.taak='FRT';
elseif taakinvoer=='p'
    mlib('write', vars.isKrachttaak, 'Data', 0);
    protocol.taak='PRT';
else
    error('Ongeldige taakinvoer')
end
    
% Maak data structure
% Hierin zit zowel het protocol als de gemeten data, zodat deze altijd
% gekoppeld blijven.
Ntrials = protocol.settrials;
totTrials = length(protocol.trials)
Nsample = settings.meetduur * settings.Fs + 1;

% Maak lege cell om geheugenruimte voor metingen te reserveren.
emptycell = cell(Ntrials, 1);
emptycell(:) = {zeros(Nsample, length(protocol.settrials))};
% Sla protocol en meetdata op in enkele structure.
dataout.protocol = protocol;
dataout.data = struct('tijd', emptycell, 'F', emptycell, 'X', emptycell);


voortgang = 0;

% Backup protocol
backupNaam = '.\autobackup\backup_protocol';
save(backupNaam, 'protocol', 'settings');
for trialKracht = 1:length(protocol.krachtniveau)
    mlib('write', vars.boodschap, 'data', 7);
    pause(2)
for trialNr = startMetTrialNr:Ntrials

    % Laad huidige trial
%     currentTrial = protocol.trials(trialNr)
    currentTrial = protocol.trials(trialNr + ((trialKracht - 1) * protocol.settrials))
    krachtniveau = protocol.krachtniveau(trialKracht);
    % Voer trial uit
    [tijd,X,F] = DoeTrial(currentTrial, vars, settings, protocol, krachtniveau);

    % Sla de gemeten data op.
    dataout.data(trialKracht).tijd(:,trialNr) = tijd;
    dataout.data(trialKracht).F(:,trialNr) = F;
    dataout.data(trialKracht).X(:,trialNr) = X;
%     dataout.data(trialNr).tijd(:) = tijd;
%     dataout.data(trialNr).F(:) = F;
%     dataout.data(trialNr).X(:) = X;

    % Sla een reservekopie van dataout op voor eventuele vastlopers of
    % problemen. Met deze data en het protocol is alle meetdata te
    % reconstrueren.
    backupData = dataout.data(trialNr)
    backupNaam = ['.\autobackup\backup_trial_', num2str(trialNr,'%05d')];
    save(backupNaam, 'backupData');
    
    % Bereken voortgang
%     voortgang = (trialNr - startMetTrialNr + 1) / (Ntrials - startMetTrialNr + 1);
    voortgang = (trialNr - startMetTrialNr + 1 + (trialKracht - 1) * protocol.settrials) / (totTrials - startMetTrialNr + 1) 
    mlib('write', vars.voortgang, 'data', voortgang)
    
    disp(['Voortgang: ', num2str(trialNr), ' / ', num2str(Ntrials)])

end
end
% Sla de meetdata op
disp('Data wordt opgeslagen')
% saveNaam = ['.\data\', datestr(now,'yyyy-mm-dd HH-MM-SS'), '_', protocol.proefpersoon, '_', num2str(protocol.stijfheid), '_', protocol.taak];
saveNaam = ['.\data\', protocol.proefpersoon, '_', protocol.taak];

save(saveNaam, 'dataout')
disp('Data is opgeslagen')

% Geef aan dat experiment klaar is
mlib('write', vars.boodschap, 'data', 4)
mlib('write',vars.Flevel,'Data',0)
mlib('write',vars.Koff,'Data',25/100)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           INIT PROPRIO                        %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitProprio
% Load
disp('Loading application to dSpace')
disp(' ')
ok=dos('scoutcmd .\controller\controller_no_k_Flevel\controller_no_k_Flevel','-echo');
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
vars.X = mlib('GetTrcVar','Model Root/offset sum X/Out1');
vars.F = mlib('GetTrcVar','Model Root/ForceFilt/Out1');
vars.stijfheid = mlib('GetTrcVar','Model Root/P:Stijfheid_setpoint.Value');
vars.Fafwijking = mlib('GetTrcVar','Model Root/Stijfheid/P:Fafwijking.Value');
vars.isKrachttaak = mlib('GetTrcVar','Model Root/Visualisatie/P:isKrachttaak.Value');
vars.visueel = mlib('GetTrcVar','Model Root/Visualisatie/P:VisSwitch.Gain');
vars.terug = mlib('GetTrcVar','Model Root/P:Terug.Value');
vars.triggerAanUit = mlib('GetTrcVar','Model Root/P:triggerAanUit.Gain');
vars.trigger = mlib('GetTrcVar','Model Root/triggerAanUit/Out1');
vars.resetRecordLight = mlib('GetTrcVar','Model Root/P:reset TriggerSensor.Value');
vars.boodschap = mlib('GetTrcVar','Model Root/P:Boodschap.Value');
vars.voortgang = mlib('GetTrcVar','Model Root/P:Voortgang.Value');
vars.Flevel = mlib('GetTrcVar','Model Root/P:Flevel.Value');
vars.Ksafe = mlib('GetTrcVar','Model Root/Safety/P:Ksafe.Value');
vars.Bsafe = mlib('GetTrcVar','Model Root/Safety/P:Bsafe.Value');
vars.Safegain = mlib('GetTrcVar','Model Root/Safety/P:safetygain.Gain');
vars.Xdoel = mlib('GetTrcVar','Model Root/P:Xdoel.Value');
vars.Koff = mlib('GetTrcVar','Model Root/Safety/P:Koff.Gain');
end

% mlib('write',vars.Safegain,'Data',1)
% mlib('write',vars.Safegain,'Data',0)
% mlib('write',vars.Flevel,'Data',0)
% mlib('write',vars.Ksafe,'Data',0)
% mlib('write',vars.Bsafe,'Data',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Offset                              %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = OffSet(vars, settings)

% Maak additionele variabelen aan
add.rawdX = mlib('GetTrcVar','Model Root/gain dX/Out1');
add.rawX = mlib('GetTrcVar','Model Root/gain X/Out1');
add.rawF = mlib('GetTrcVar','Model Root/gain Fc/Out1');
add.offsetdX = mlib('GetTrcVar','Model Root/P:offset dX.Value');
add.offsetX = mlib('GetTrcVar','Model Root/P:offset X.Value');
add.offsetF = mlib('GetTrcVar','Model Root/P:offset Fc.Value');
add.demping = mlib('GetTrcVar','Model Root/P:Damping.Value');

% Zet huidige offsets op nul
mlib('write', add.offsetX, 'Data', 0);
mlib('write', add.offsetdX, 'Data', 0);
mlib('write', add.offsetF, 'Data', 0);

% Lichte demping en stijfheid om naar absolute nulpunt te komen
% mlib('write', vars.stijfheid, 'Data', 50);
% mlib('write', vars.Fafwijking, 'Data', 0);
% mlib('write', add.demping, 'Data', 5);

mlib('write', vars.Koff, 'Data', 50/100);

input('Check of de hardwarematige massa, stijfheid en demping op nul staan en druk enter','s')
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
    'TraceVars',[add.rawX; add.rawdX; add.rawF],...
    'AcquisitionMode','single_shot',...
    'DownSampling',downSample,...
    'NumSamples',nrDatSamples);

% start capture (position only) and wait till ready
mlib('StartCapture');

while (mlib('CaptureState')~=0)
    pause(0.01)
end

output = mlib('FetchData');
offsetX = mean(output(1,:)) + settings.nulpositie;
offsetX = sign(offsetX) * min(abs(offsetX), 0.2);
offsetdX = mean(output(2,:));

% Stel positie offset waarden in op dspace
mlib('write', add.offsetX, 'Data', offsetX);
mlib('write', add.offsetdX, 'Data', offsetdX);

% Wacht totdat nieuwe positie gevonden is
pause(4)

% start capture (force) and wait till ready
mlib('StartCapture');

while (mlib('CaptureState')~=0)
    pause(0.01)
end

output = mlib('FetchData');
offsetF = mean(output(3,:));

% Stel offset waarden in op dspace
mlib('write', add.offsetF, 'Data', offsetF);
mlib('write',vars.Safegain,'Data',1)
mlib('write', vars.Koff, 'Data', 0);
disp('Offsetten: klaar.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Doe enkele trial meting             %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tijd,X,F] = DoeTrial(trial, vars, settings, protocol, krachtniveau)

% Duur van meting instellen
nrDatSamples = settings.meetduur * settings.Fs + 1;
downSample = settings.Fds / settings.Fs;

% Initialize MTRACE for data acquisition
mlib('Set',...
    'Trigger','On',...
    'TriggerVariable', vars.trigger,...
    'TriggerLevel', 0.1,...
    'TriggerEdge', 'rising',...
    'TraceVars',[vars.X; vars.F],...
    'AcquisitionMode','single_shot',...
    'DownSampling',downSample,...
    'NumSamples',nrDatSamples);

% Wacht tot X in beginpositie
disp('Beweeg handle naar beginpositie')
mlib('write', vars.boodschap, 'data', 1)
mlib('write', vars.terug, 'data', 1)
Xnow = mlib('read', vars.X);
while Xnow > 0.002
    pause(0.1)
    Xnow = mlib('read', vars.X);
end
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

mlib('write', vars.terug, 'data', 0)

% Stel trialeigenschappen in
% mlib('write', vars.stijfheid, 'data', protocol.stijfheid);
% mlib('write', vars.Fafwijking, 'data', trial.Fafwijking);
mlib('write', vars.visueel, 'data', trial.visueel);
% mlib('write', vars.Flevel, 'data', 1);
% pause(1)
mlib('write', vars.Flevel, 'data', krachtniveau);

if strcmp(protocol.taak,'kracht')
    if trial.visueel
        mlib('write', vars.boodschap, 'data', 2);
    else
        mlib('write', vars.boodschap, 'data', 3);
    end
elseif strcmp(protocol.taak,'positie')
    if trial.visueel
        mlib('write', vars.boodschap, 'data', 5);
    else
        mlib('write', vars.boodschap, 'data', 6);
    end
else
    error('taak onduidelijk')
end

% Zet trigger op scherp
mlib('write', vars.triggerAanUit, 'data', 1)

% Start getriggerde meting
disp('Meting trial gestart.')
mlib('StartCapture');
while (mlib('CaptureState')~=0)
    pause(0.1)
end
pause(0.6)
% Reset triggersensor (tbv record-lampje)
mlib('write', vars.resetRecordLight, 'data', 1)
pause(0.1)
mlib('write', vars.resetRecordLight, 'data', 0)

% Zet trigger uit en geef aan dat proefpersoon terug moet
mlib('write', vars.triggerAanUit, 'data', 0)
mlib('write', vars.terug, 'data', 1)

% Lees data uit
output = mlib('FetchData');
disp('Meting trial klaar.')
X = output(1,:)';
F = output(2,:)';
tijd = (0:nrDatSamples-1)' / settings.Fs;
end
