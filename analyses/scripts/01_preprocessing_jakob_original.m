%%
addpath('/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0')
addpath('/home/marie/Documents/MATLAB/NoiseTools')
addpath('/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0/plugins/xdfimport1.18/xdf')
eeglab;
cd('/home/marie/Documents/uni/8_SS_21/mbensienBA/data/02_main');

%% check for noisetools
try, nt_greetings; catch, disp('You must download NoiseNools from http://audition.ens.fr/adc/NoiseTools/'); return; end
%%
%loading xdf into mobilab 
allDataStreams = dataSourceXDF('pilot2_081220.xdf','MobiLabData');
EEG = allDataStreams.item{9};

%% export into eeglab (8 is the eeg channel index) and save as set and name the set
%addpath('./eeglab2020_0');
EEG = allDataStreams.export2eeglab( 8 )
EEG = pop_saveset( EEG, 'filename','newest_pilot2','filepath','E:\\Uni\\Bachelor\\');
EEG = pop_editset(EEG, 'setname', 'newest_pilot2', 'run', []);
EEG = eeg_checkset(EEG)
%% load the before mentioned .set into eeglab
sub = 'pilot1';
filename = 'new_pilot1.set';
filepath = 'E:/Uni/Bachelor';
EEG = pop_loadset(filename, filepath);
%%
EEG.setname = filename;
setname = EEG.setname;
EEG = eeg_checkset(EEG);
%% STEP 2: FILTER, RESAMPLE

% List of Variables
srate = EEG.srate/4; %from 1024 to 256 srate
HP_cutoff = 0.1;
LP_cutoff = 120;

EEG.preprocessing = [];

% Downsampling to 256 Hz
EEG = pop_resample(EEG, srate);
EEG.preprocessing = [EEG.preprocessing 'Resampled,'];
EEG.comments = pop_comments(EEG.comments,'','Dataset was downsampled to 256 Hz.',1);

% Low pass filter under 120 Hz
EEG = pop_eegfiltnew(EEG, [], LP_cutoff);
EEG.preprocessing = [EEG.preprocessing 'Lowpass,'];
EEG.comments = pop_comments(EEG.comments,'','Dataset was lowpass filtered at 120 Hz.',1);

% High pass filter over 0.1Hz
EEG = pop_eegfiltnew(EEG, HP_cutoff, []);
EEG.preprocessing = [EEG.preprocessing 'Highpass,'];
EEG.comments = pop_comments(EEG.comments,'','Dataset was highpass filtered at 0.1 Hz.',1);

%% Zapline filter out 50 Hz line noise
%parameters 
FLINE = 50/256; % line frequency
nremove = 2;
d = permute(EEG.data, [2,1]);
d = nt_zapline(d ,FLINE, nremove);
EEG.data = permute(d, [2,1]);
EEG.preprocessing = [EEG.preprocessing 'Zapline,'];
EEG.comments = pop_comments(EEG.comments,'','Dataset was zapline filtered at 50 Hz.',1);
%%
pop_eegplot(EEG);
%% Locations
% TODO: Import channel locations
%EEG = pop_select( EEG, 'nochannel',{'BIP65','BIP66','BIP67','BIP68','AUX69','AUX70', 'AUX71', 'AUX72'});
eeglabpath = 'E:/Uni/Bachelor/eeglab2020_0';
EEG=pop_chanedit(EEG, 'lookup','E:\\Uni\\Bachelor\\eeglab2020_0\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','load',{'E:\\Uni\\Bachelor\\eeglab2020_0\\sample_locs\\my_locs.ced','filetype','autodetect'});
pop_spectopo(EEG);
%%
% save this dataset
% change the setname from pilot01 to pilot01_resampled_filtered
EEG = pop_editset(EEG, 'setname', sprintf('newest_pilot1_zapline2test',setname));
EEG = pop_saveset(EEG, 'filename',sprintf('newest_pilot1_zapline2test',setname),'filepath',filepath);

%% Reject parts of continous data and save them
winlength = 10;                 % time range to display in seconds
scale = 70;                     % amplitude range for channels (uV)



    eegplot(EEG.data,'command','rej=TMPREJ','title',['Reject time intervals manually'],'events',EEG.event,'spacing',scale,'eloc_file',EEG.chanlocs);
%% save information about removed parts of data
    tmprej = eegplot2event(rej, -1);
%%  Save file and rejected data parts  
filepath = 'E:/Uni/Bachelor';
    EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
    EEG.preprocessing = [EEG.preprocessing 'Cleaning,'];
    EEG = pop_editset(EEG, 'setname', sprintf('newest_pilot1_Clean_Zapline2',setname));
    EEG = pop_saveset(EEG, 'filename',sprintf('newest_pilot1_Clean_Zapline2',setname),'filepath', filepath);
                
    save(fullfile(filepath,sprintf('%s_cleaningTimes.mat',setname)),'tmprej','rej');
%% save set
EEG = pop_editset(EEG, 'setname', sprintf('newest_pilot1_Clean',setname));
EEG = pop_saveset(EEG, 'filename',sprintf('newest_pilot1_Clean',setname),'filepath', filepath);
save(fullfile(filepath,sprintf('newest_pilot1_cleaningTimes.mat',setname)),'tmprej','rej');
%% ICA    
addpath('./amica')
outDir = fullfile(filepath,'amica');
%for the ICA computation, highpass-filter data at 2Hz to get better
%Eye-Components
EEG2Hz = pop_eegfiltnew(EEG, 2, []);   % highpass
runamica15(EEG2Hz.data,'outdir',outDir);

%% Apply ICA weights to data
%load ICA results
addpath([fullfile(filepath) '/amica'])
                icapath = fullfile(filepath,'amica');
                mod = loadmodout15(icapath);
                

                %%% APPLY ICA WEIGHTS TO DATA %%%
                EEG.icasphere = mod.S;
                EEG.icaweights = mod.W;
                EEG = eeg_checkset(EEG);
                EEG.preprocessing = [EEG.preprocessing 'AMICA,'];
%% 
EEG = pop_editset(EEG, 'setname', 'newest_pilot1_AMICA_proper_Zapline2', 'run', []);
EEG = pop_saveset( EEG, 'filename','newest_pilot1_AMICA_proper_Zapline2','filepath','E:\\Uni\\Bachelor\\');
EEG = eeg_checkset(EEG);
%%
pop_eventstat(EEG)
%% plot components
% first compute ica weights
EEG.icaact = eeg_getdatact(EEG,'component', 1:size(EEG.icaweights, 1));
%% then plot the activation of the components
eegplot(EEG.icaact([8],:,:), 'srate', EEG.srate)