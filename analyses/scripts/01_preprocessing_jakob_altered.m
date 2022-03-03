%%% This script was altered from Jakob's original script and was only used
%%% for testing. 
%% Add all paths to dependencies: eeglab, NoiseTools, xdfimport
% start eeglab and go to data path
addpath('/home/marie/Documents/MATLAB/eeglab2020_0')  % make sure to use 2020 version 
addpath('/home/marie/Documents/MATLAB/NoiseTools')
addpath('/home/marie/Documents/MATLAB/eeglab2020_0/plugins/xdfimport1.18/xdf')
eeglab;
cd('/home/marie/Documents/uni/8_SS_21/mbensienBA/data/0_raw');

%% check for noisetools
try, nt_greetings; catch, disp('You must download NoiseNools from http://audition.ens.fr/adc/NoiseTools/'); return; end

%% Load data into MobiLab, export EEG streams into eeglab set, save 
allDataStreams = dataSourceXDF('/home/marie/Documents/uni/8_SS_21/mbensienBA/data/0_raw/10_r3_140621.xdf','MobiLabData'); % replace 04_r1_020621 with current filename
% 9 is the eeg stream index in xdf data frame
EEG = allDataStreams.export2eeglab(9);
EEG = pop_saveset( EEG, 'filename','10_r3_140621','filepath','/home/marie/Documents/uni/8_SS_21/mbensienBA/data/1_preprocessed');
EEG = pop_editset(EEG, 'setname', '10_r3_140621', 'run', []);
EEG = eeg_checkset(EEG); % checks consistency of set
filename = '10_r3_140621.set';
filepath = '/home/marie/Documents/uni/8_SS_21/mbensienBA/data/1_preprocessed';
EEG = pop_loadset(filename, filepath);
EEG.setname = filename;
setname = EEG.setname;
EEG = eeg_checkset(EEG); % checks consistency of set


%% STEP 1: FILTER, RESAMPLE

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

%% Load set if this is where you set off
filename = '05_r3_080621_eeg_filtered.set';
filepath = '/home/marie/Documents/uni/8_SS_21/mbensienBA/data/00_raw/rejected';
EEG = pop_loadset(filename, filepath);
EEG.setname = filename;
setname = EEG.setname;
EEG = eeg_checkset(EEG); % checks consistency of set

%% Import electrode locations
% Exclude additional channels which were not used and import channel locations
EEG = pop_select(EEG, 'nochannel',{'BIP65','BIP66','BIP67','BIP68','AUX69','AUX70', 'AUX71', 'AUX72'});
eeglabpath = '/home/marie/Documents/MATLAB/eeglab2020.0';
EEG = pop_chanedit(EEG, 'lookup','/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0/plugins/dipfit4.0/standard_BESA/standard-10-5-cap385.elp','load',{'/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0/sample_locs/my_locs.ced','filetype','autodetect'})
% save set with channel labels and locations
EEG.preprocessing = [EEG.preprocessing 'chanlocs, empty channels excluded,'];
EEG.comments = pop_comments(EEG.comments,'','Channel locations and labels added. Empty channels removed.',1);
EEG = pop_saveset( EEG, 'filename','10_r3_filtered_chanlocs.set','filepath','/home/marie/Documents/uni/8_SS_21/mbensienBA/data/1_preprocessed');

%% Inspect data and remove noisy channels
% Look at power spectrum across electrodes and topographical distribution at key frequencies
pop_spectopo(EEG); % sample 100 percent, frequencies with peaks in range 2-50

% Find names of electrodes by index 
electrodeNums = [13, 28, 62, 33];
electrodeNames = strings(size(electrodeNums));
for electrode=1:length(electrodeNums)
    electrodeNames(electrode) = EEG.chanlocs(electrodeNums(electrode)).labels;
end
electrodeNames
% or look up a single electrode name by index
EEG.chanlocs(13).labels 

% Look at continuous data and reject noisy channels
pop_eegplot(EEG); 
EEG = pop_select(EEG, 'nochannel', {'M2'}); 
EEG.preprocessing = [EEG.preprocessing 'M2 rejected,'];

%% Save changes to set 
EEG = pop_editset(EEG, 'setname', sprintf('05_r3_filtered_chanlocs_exclM2', setname));
EEG = pop_saveset(EEG, 'filename',sprintf('05_r3_filtered_chanlocs_exclM2', setname),'filepath',filepath);


%% Re-reference to Cz
EEG = pop_reref( EEG, 16);
EEG = pop_newset(ALLEEG, EEG, 1,'setname','09_r1_110621_eeg_filtered_channelReject_referenceCz','savenew','09_r1_110621_eeg_filtered_channelReject_referenceCz','comments',strvcat('Dataset was downsampled to 256 Hz.','Dataset was lowpass filtered at 120 Hz.','Dataset was highpass filtered at 0.1 Hz.','Dataset was zapline filtered at 50 Hz.','Dataset was referenced to Cz. Cz channel stays in data.'),'gui','off'); 
EEG = eeg_checkset( EEG );

%% Manually reject parts of continous data and save them
winlength = 10;                 % time range to display in seconds
scale = 70;                     % amplitude range for channels (uV)

eegplot(EEG.data,'command','rej=TMPREJ','title',['Reject time intervals manually'],'events',EEG.event,'spacing',scale,'eloc_file',EEG.chanlocs);
% save information about removed parts of data
    tmprej = eegplot2event(rej, -1);
%  Save file and rejected data parts  
filepath = '/home/marie/Documents/uni/8_SS_21/mbensienBA/data/01_preprocessed';
    EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
    EEG.preprocessing = [EEG.preprocessing 'Cleaning,'];
    EEG = pop_editset(EEG, 'setname', sprintf('04_r1_020621_Clean_Zapline2',setname));
    EEG = pop_saveset(EEG, 'filename',sprintf('04_r1_020621_Clean_Zapline2',setname),'filepath', filepath);
                
    save(fullfile(filepath,sprintf('%s_cleaningTimes.mat',setname)),'tmprej','rej');

    %% Save set
EEG = pop_editset(EEG, 'setname', sprintf('04_r1_020621_Clean',setname));
EEG = pop_saveset(EEG, 'filename',sprintf('04_r1_020621_Clean',setname),'filepath', filepath);
save(fullfile(filepath,sprintf('04_r1_020621_cleaningTimes.mat',setname)),'tmprej','rej');
EEG = pop_loadset('04_r1_020621_Clean.set', filepath);

%% Independent component analysis
% high-pass filter at 2Hz for better eye components
addpath('/home/marie/Documents/MATLAB/eeglab2020_0/plugins/AMICA1.5.2')
outDir = fullfile(filepath,'01_amica');
EEG2Hz = pop_eegfiltnew(EEG, 2, []);   % highpass
runamica15(EEG2Hz.data,'outdir',outDir);


%% Apply ICA weights to data
%load ICA results
addpath([fullfile(filepath) '/01_amica'])
icapath = fullfile(filepath,'01_amica');
mod = loadmodout15(icapath);
                

%%% APPLY ICA WEIGHTS TO DATA 
EEG.icasphere = mod.S;
EEG.icaweights = mod.W;
EEG = eeg_checkset(EEG); % checks consistency of set
EEG.preprocessing = [EEG.preprocessing 'AMICA,'];

%% Save changes to set 
% after ICA
EEG = pop_editset(EEG, 'setname', '04_r1_020621_AMICA_proper_Zapline2', 'run', []);
EEG = pop_saveset( EEG, 'filename','04_r1_020621_AMICA_proper_Zapline2','filepath','/home/marie/Documents/uni/8_SS_21/mbensienBA/data/preprocessed/02_ICA');
EEG = eeg_checkset(EEG); % checks consistency of set

%% look this up?
pop_eventstat(EEG)

%% plot components
% first compute ica weights
EEG.icaact = eeg_getdatact(EEG,'component', 1:size(EEG.icaweights, 1));
% then plot the activation of the components
eegplot(EEG.icaact([8],:,:), 'srate', EEG.srate)