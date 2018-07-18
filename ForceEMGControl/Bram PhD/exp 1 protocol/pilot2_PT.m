%% Dit is een voorbeeldexperiment.
clear all
clc

%% Krachtniveau instellen

% Krachtniveau = input('Krachtniveau: ');

% Krachtlevels = [10:20:170];
Krachtlevels = [10 20 30 40 50 60 70 80];
randomkracht = randperm(length(Krachtlevels));
Krachtniveau = Krachtlevels(randomkracht); 

%% Settings
settings.meetduur = 0.3;
settings.nulpositie = -0.150;
settings.Fds=5000;	% dspace op 5000 Hz (simulink model sample fq = 500Hz)
settings.Fs=500;   	% sampling op 100 Hz

%% Definieer twee standaard trials en maak daaruit een serie.

% Definieer Normale trial
trNorm.omschrijving = 'norm';
%trNorm.stijfheid = Stijfheid;
trNorm.Fafwijking = 0;
trNorm.visueel = 0;

% Definieer train trial
trTrain.omschrijving = 'train';
%trTrain.stijfheid = Stijfheid;
trTrain.Fafwijking = 0;
trTrain.visueel = 1;

% Maak de serie met trials aan
% Eerst 15 training trials, dan om de andere trial een normale
clear alleTrials
alleTrials = [];
krachtTrials (1:21,1) = trTrain;
krachtTrials (6:2:end) = trNorm; 

for i = 1:length(Krachtlevels)
    alleTrials = [alleTrials; krachtTrials]
end

% Display het resultaat, ter controle
for ii = 1:length(alleTrials);
    disp(['Trial ' num2str(ii), ': ', alleTrials(ii).omschrijving]);
end

%% Maak het protocol aan
protocol.datum = '';            % wordt automatisch ingevuld
protocol.proefpersoon = '';     % wordt naar gevraagd bij starten experiment
protocol.stijfheid = '0';
protocol.krachtniveau = Krachtniveau;
protocol.omschrijving = 'number of repetitions';
protocol.trials = alleTrials;
protocol.settrials = length(krachtTrials);

% Display het resultaat, ter controle
protocol



%% Start het experiment
startMetTrialNr = 1; % Maakt het mogelijk een afgebroken protocol te vervolgen.
runprotocol(protocol, settings, startMetTrialNr)
