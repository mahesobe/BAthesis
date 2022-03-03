% This script allows you to preprocess your EEG data
% semiautomatically.
% It contains the following steps:
% 1 - loading, filtering, downsampling
% 2 - channel position, reference to Cz, channel cleaning
% 3 - cleaning continous data
% 4 - temporarily apply 2 Hz high-pass filter and calculate ICA with AMICA
% 5 - apply AMICA weights 
% 6 - show AMICA components, clean them (from continous data)
% 7 - interpolation of missing channels
%
%
% The script developed over time, but mostly it was written by @tkietzma (automation, GUI)
% @sitimm (adjusting, behavioral parts), @agert (adjusting, making it ready for EEG training)
% + in this case Marie Bensien who used it for her Bachelor's thesis
%% to restart everything
% do not run when running it automatically
 clear variables
 close all;
 clc;
 restoredefaultpath
%% Load necessary code
cd /home/marie/Documents/uni/8_SS_21/mbensienBA/data/ %folder to save matlab steps
savepath = "/home/marie/Documents/uni/8_SS_21/mbensienBA/data";

% add all necesarry tools
addpath('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/preprocessing_helpers');
addpath('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/NoiseTools');
remove_from_path('eeglab');  % remove old eeglab
addpath('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/eeglab2020_0');
addpath('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/eeglab2020_0/plugins/xdfimport1.18/xdf');

eeglab;

% Where is your data stored?
basepath='/home/marie/Documents/uni/8_SS_21/mbensienBA/data/01_raw';

autorun = 0;
filtType   ='acausal';	% how do you want to filter? Options: acausal or causal
filtHigh   = 0.1;       % which frequency do you want to use as high pass for filtering? Recomended: 0.1Hz
filtLow    = 120;       % which frequeny do you want to use as low pass for filtering?
eyetracking= 0;         % did you eyetracking and want to integrate that data? This is not implemented in the EEG_Training but IntoTheWild
resampleFreq = 512;     

%%% LOAD RAW DATA %%%
%select data file to work on
cntpath = fullfile(basepath,'10_r3_140621.xdf');
if exist(cntpath,'file')
    [filepath,filename,ext] = fileparts(cntpath); % define variables
    % to keep the same as the manual selection procedurezap 
    filepath = [filepath filesep];
    filename = [filename ext];
else
    cntpath = fullfile(basepath,sprintf('VP%u',sub),'*.*');
    [filename, filepath]=uigetfile(cntpath);
end

setname=dir(fullfile(filepath,filename));
setname=setname.name(1:end-4); % to get the name without the extension

%% Set Up GUI

try
    %temp:
    sub = setname;
    
    h=GUI; %script in 'preprocessing_helpers'
    data = guihandles(h);
    set(data.titletext, 'String', sprintf('Subject %u, Task %s',sub,setname)); % add name to gui
    
    
    test_preproc_finished;
    if autorun == 0 % manual mode
        uiwait(h)
        selection = getappdata(h, 'selection');
    end
    
    exitrequest=false;
    while ~exitrequest
        
        
        if autorun ~= 0
            selection = auto_step_list(autorun);
            autorun = autorun + 1;
        end
        
        
        switch selection
            
            case 0
                exitrequest=1;
                close(h);
                break;
                
            case 1 %filtering
                             
                
                % create a target folder
                mkdir(fullfile(savepath,'02_preprocessed/'))
                
                % load data into Mobilab 
                allDataStreams = dataSourceXDF(fullfile(basepath, filename),'MobiLabData');
                EEG = allDataStreams.item{9};
                EEG = allDataStreams.export2eeglab( 9 ) % 9 is the stream which contains EEG signal (in rare cases, it's a different one)
                
                path = fullfile(savepath,'02_preprocessed/');
                EEG = pop_saveset( EEG, 'filename', sprintf('0_%s_raw',setname), 'filepath', char(fullfile(savepath,'02_preprocessed/')));

                
                EEG = pop_editset(EEG, 'setname', sprintf('0_%s_raw',setname));
                EEG = eeg_checkset(EEG)

                
                EEG.preprocessing = 'Raw,'; %every step gets stored in .preprocessing so you can see what happened to the data set
                                
                % make EEG file that is untouched for later use
                trEEG=EEG;
                
                
                %%% HighPass FILTER %%%
                % filter out low frequencies
                if strcmp(filtType,'acausal')
                    EEG = pop_eegfiltnew(EEG, filtHigh, []);
                elseif strcmp(filtType,'causal')
                    EEG = pop_eegfiltnew(EEG, filtHigh,[],[],[],[],0,1);
                end
                EEG.preprocessing = [EEG.preprocessing 'Highpass' num2str(filtHigh) ', '];
                
                %%% LowPass FILTER %%%
                % filter out high frequencies
                EEG = pop_eegfiltnew(EEG, [], filtLow);
                EEG.preprocessing = [EEG.preprocessing 'Lowpass ' num2str(filtLow) ', '];
                
                %%% Zapline FILTER %%%
                % filter out 50 Hz line noise 
                FLINE = 50/512; % line frequency
                nremove = 2;
                d = permute(EEG.data, [2,1]);   
                d = nt_zapline(d ,FLINE, nremove); 
                EEG.data = permute(d, [2,1]);
                EEG.preprocessing = [EEG.preprocessing 'Zapline 50, '];
                    
                %%% DOWNSAMPLING %%%
                % lower sampling rate from default (1024Hz) to 512Hz
                EEG = pop_resample(EEG, resampleFreq);
                EEG.preprocessing = [EEG.preprocessing 'Resample ' num2str(resampleFreq) ', '];
                trEEG = pop_resample(trEEG, resampleFreq);
                
                
                % Save the data
                EEG = pop_editset(EEG, 'setname', sprintf('1_%s_filtered_resample',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('1_%s_filtered_resample',setname),'filepath',char(fullfile(savepath,'02_preprocessed/')));
                
                trEEG = pop_saveset(trEEG, 'filename',sprintf('%s_triggers',setname),'filepath',char(fullfile(savepath,'02_preprocessed/')));
                
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')

                
            case 2 %behav, position, channelselect
                
                fprintf('Loading filtered and resampled files...\n')
                
                %load files
                EEG = pop_loadset(sprintf('1_%s_filtered_resample.set',setname),char(fullfile(savepath,'02_preprocessed/')));
                trEEG = pop_loadset(sprintf('%s_triggers.set',setname),char(fullfile(savepath,'02_preprocessed/')));
                
                
                %% Channel positions
                % drop empty channels, add locations and reference to Cz 
                EEG = pop_select(EEG, 'nochannel',{'BIP65','BIP66','BIP67','BIP68','AUX69','AUX70', 'AUX71', 'AUX72'});
                eeglabpath = '/home/marie/Documents/MATLAB/eeglab2020.0';
                EEG = pop_chanedit(EEG, 'lookup','/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0/plugins/dipfit4.0/standard_BESA/standard-10-5-cap385.elp','load',{'/home/marie/Documents/MATLAB/eeglab_current/eeglab2021.0/sample_locs/my_locs.ced','filetype','autodetect'})
                EEG.preprocessing = [EEG.preprocessing 'ChanPos,'];
                EEG = pop_reref(EEG, 16);
                EEG.preprocessing = [EEG.preprocessing 'RefCz,'];
                EEG_tmp=EEG;                
                
                EEG = pop_editset(EEG, 'setname', sprintf('2_%s_chanlocsRefCZ',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('2_%s_chanlocsRefCZ',setname),'filepath',fullfile(char(fullfile(savepath,'02_preprocessed/'))));
                

                %%% CHANNEL CLEANING %%%
                % which channels do not contain EEG data by default?
                alldel = {'BIP2' 'BIP3' 'BIP4' 'BIP5' 'BIP6' 'BIP7' 'BIP8' 'AUX1' 'AUX2' 'AUX3' 'AUX4' 'AUX5' 'AUX6' 'AUX7' 'IZ' 'TIME' 'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'R-GAZE-X' 'R-GAZE-Y' 'R-AREA' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA', 'R_GAZE_X' 'R_GAZE_Y' 'R_AREA' 'INPUT'};
                
                str=[];
                delindex =[];
                for k=1:numel(EEG.chanlocs)
                    str{k}=EEG.chanlocs(k).labels;
                    
                    if strmatch(str{k},alldel)
                        delindex(k)=1;
                    else
                        delindex(k)=0; 
                    end
                end
                delindex = find(delindex);
                targetchannels = str(setdiff(1:length(str),delindex));
                
                
                
                for k=1:numel(EEG.chanlocs)
                    str{k}=EEG.chanlocs(k).labels;   
                end
                
                % If you already cleaned the channels, recover that information
                if exist(fullfile(char(fullfile(savepath,'02_preprocessed/')),sprintf('2_%s_channelrejTriggersXensor.mat',setname)))

                    tmp=load(fullfile(savepath,'02_preprocessed',sprintf('2_%s_channelrejTriggersXensor.mat',setname)));
                    
                    %make sure unique channels are selected
                    tmp.chan_del=unique(tmp.chan_del);
                    
                    
                    fprintf('PREVIOUS CHANNEL CLEANING FOUND\n')
                    fprintf('MARKED CHANNELS ARE: ')
                    tmp.chanids=[];
                    
                    for b=[1:length(tmp.chan_del)]
                        fprintf('%s,',char(tmp.chan_del(b)))
                        if ~isempty(strmatch(tmp.chan_del(b),str,'exact'))
                            tmp.chanids =[tmp.chanids strmatch(tmp.chan_del(b),str,'exact')];
                        end
                    end
                    fprintf('\n\n');
                    delindex = unique([delindex tmp.chanids]);
                end
                
                
                % automatically remove the default channels from the data set that will be visualized in the GUI
                tmpEEG = pop_select(EEG, 'nochannel', alldel);

                if autorun == 0 % manual mode
                    eegplot(tmpEEG.data,'srate',tmpEEG.srate,'eloc_file',tmpEEG.chanlocs,'events',tmpEEG.event)
                    
                    %wait for user to note channels
                    choice = menu('Noisy channels noted?','yes');
                    
                    
                    %create menu for input
                    [s,v] = listdlg('PromptString','Select Channels:',...
                        'SelectionMode','multiple',...
                        'InitialValue',delindex,...
                        'ListString',str);
                    s=[s delindex];
                    
                else
                    s = [delindex]; % only the previously loaded ones
                end
                
                chan_del = str(s);
                
                %remove channels from EEG set
                EEG = pop_select(EEG, 'nochannel', chan_del);
                
                %keep full chanlocs for later reinterpolation
                complete_chanlocs = EEG_tmp.chanlocs;
                
                EEG.preprocessing = [EEG.preprocessing 'ChannelReject,'];
                
                % file
                EEG = pop_editset(EEG, 'setname', sprintf('2_%s_channelrejTriggersXensor',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('2_%s_channelrejTriggersXensor',setname),'filepath',fullfile(char(fullfile(savepath,'02_preprocessed/'))));
                
                save(fullfile(char(fullfile(savepath,'02_preprocessed/')),sprintf('2_%s_channelrejTriggersXensor.mat',setname)),'chan_del', 'targetchannels','complete_chanlocs');
                % look at data again to confirm 
                eegplot(EEG.data,'srate',EEG.srate,'eloc_file', EEG.chanlocs,'events',EEG.event)


                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
                
                
            case 3 %data cleaning
                fprintf('Loading channel cleaned files...\n')
                
                
                %load files
                EEG = pop_loadset(sprintf('2_%s_channelrejTriggersXensor.set',setname), fullfile(char(fullfile(savepath,'02_preprocessed/'))));
                if isempty(EEG.urevent) %creates urevent structure for later use
                    EEG = eeg_checkset(EEG,'makeur');
                end
                
                
                %to exclude luminance sensor channels from GUI (data not lost, but excluded from analysis)
                auxdel = find(cellfun(@(x)~isempty(x),strfind({EEG.chanlocs(:).labels},'AUX')));
                tmpcleanEEG= pop_select(EEG, 'nochannel', auxdel);
                
                
                %%% CLEANING %%%
                %using data scroll, mark all regions that need to be excluded. Once
                %finished, click "reject". Make sure the cleaning times from 'eegh' match the stored ones shown by running the code below.
                rej = [];
                if exist(fullfile(filepath,'02_preprocessed',sprintf('3_%s_cleaningTimes.mat',setname)),'file')
                    load(fullfile(filepath,'02_preprocessed',sprintf('3_%s_cleaningTimes.mat',setname)),'tmprej','rej');
                    
                    %make sure the number of channels in rej matches the data
                    rej = [rej(:,1:5) zeros(size(rej,1),tmpcleanEEG.nbchan)];
                    
                    if autorun == 0 % manual mode
                        eegplot(tmpcleanEEG.data,'winrej',rej,'command','rej=TMPREJ;','srate',tmpcleanEEG.srate,'eloc_file',tmpcleanEEG.chanlocs,'events',tmpcleanEEG.event);
                    end
                else
                    eegplot(tmpcleanEEG.data,'command','rej=TMPREJ;','srate',tmpcleanEEG.srate,'eloc_file',tmpcleanEEG.chanlocs,'events',tmpcleanEEG.event);
                end
                
                if autorun == 0 % manual mode
                    choice=menu('Cleaning finished?','yes');
                end
                
                
                %transform rejected time windows to EEGlab interpretable structure
                tmprej = eegplot2event(rej, -1);
                
                
                %reject marked parts
                EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
                EEG.preprocessing = [EEG.preprocessing 'Cleaning,'];
                
                
                %save files
                EEG = pop_editset(EEG, 'setname', sprintf('3_%s_Clean',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('3_%s_Clean',setname),'filepath',fullfile(char(savepath),'02_preprocessed/'));
                
                save(fullfile(char(savepath),'02_preprocessed/',sprintf('3_%s_cleaningTimes.mat',setname)),'tmprej','rej');
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 4 %ICA
                fprintf('Loading data cleaned files...\n')
                
                %make sure that you don't rerun ICA automatically
                if autorun ~= 0 % manual mode
                    error('you should not have reached here. No new ICA should be started during autorun-mode')
                end
                
                if exist(fullfile(char(savepath),'02_preprocessed/amica'))
                    fprintf('ICA folder found')
                    %wait for to confirm that ICA shall be rerun
                    choice = menu('Do you really want to rerun the ICA?','yes');
                end
                
                
                %load files
                EEG = pop_loadset(sprintf('3_%s_Clean.set',setname),fullfile(char(savepath),'02_preprocessed/'));
                
                
                %to exclude luminance sensor channels from ICA (data not lost, but excluded from analysis)
                auxdel = find(cellfun(@(x)~isempty(x),strfind({EEG.chanlocs(:).labels},'AUX')));
                EEG = pop_select(EEG, 'nochannel', auxdel);
                
                
                %for the ICA computation, highpass-filter data at 2Hz to get better
                %Eye-Components
                EEG2Hz = pop_eegfiltnew(EEG, 2, []);   % highpass
                
                
                %%% ICA %%%
                addpath(fullfile('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/amica')) 
                mkdir(fullfile(char(savepath),'02_preprocessed/'),'amica')
                permission_cleanup(filepath);
                outDir = fullfile(char(savepath),'02_preprocessed/amica');
                
                
                %run AMICA locally with default parameters
                runamica12(EEG2Hz.data,'outdir',outDir,'use_queue',0,'qsubname',['ica_VP' num2str(sub)]);
                
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 5 %Apply AMICA weights 
                fprintf('Loading data cleaned files and ICA information...\n')
                
                addpath(fullfile(char(savepath),'/02_preprocessed/', '/amica'))
                addpath(fullfile('/home/marie/Documents/uni/8_SS_21/mbensienBA/analyses/01_preprocessing/Matlab-resources/amica')) 

                
                
                %load files
                EEG = pop_loadset(sprintf('3_%s_Clean.set',setname), fullfile(char(savepath),'02_preprocessed/'));
                load(fullfile(char(savepath),'02_preprocessed/',sprintf('2_%s_channelrejTriggersXensor.mat',setname)))
                
                
                %load ICA results
                icapath = fullfile(char(savepath),'02_preprocessed/','amica');
                mod = loadmodout12(icapath);
                
                
                %%% APPLY ICA WEIGHTS TO DATA %%%
                EEG.icasphere = mod.S;
                EEG.icaweights = mod.W;
                EEG = eeg_checkset(EEG);
                EEG.preprocessing = [EEG.preprocessing 'AMICA,'];
                
                
                %%% EPOCH DATA %%%
                % this epoching is done only for the vizualization later on
                % we want to keep the actual data set continous so that we are not restricted by this step later on
                %EDIT: removed bc no epochs available yet 
                %epoch_trigg={'1','2','3','4','5'};
                %window     =[-0.1 0.3];
                %baseline   =[-100 0];
                
                %[EEG,indices] = pop_epoch( EEG, epoch_trigg, window, 'epochinfo', 'yes');
                %EEG.orig_indices = indices;
                
                %EEG = pop_rmbase( EEG, baseline);
                %EEG.preprocessing = [EEG.preprocessing 'Epoched,'];
                
                
                %save files
                EEG = pop_editset(EEG, 'setname', sprintf('4_%s_ICA',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('4_%s_ICA',setname),'filepath',fullfile(char(fullfile(savepath,'02_preprocessed/'))));
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 6 %component cleaning
                %% Before this step, the dataset was loaded into the eeglab GUI,
                %% the AMICA weights applied and the components were labeled with
                %% IClabel, which allows to inspect extended properties of
                %% components. Note components you want to reject. Move on
                comps_to_rej = [];
                
                %load files
                EEG = pop_loadset(sprintf('4_%s_ICA.set',setname),fullfile(char(savepath),'02_preprocessed/'));
                EEGcontinuous= pop_loadset(sprintf('3_%s_Clean.set',setname),fullfile(char(savepath),'02_preprocessed/'));
                
                icapath = fullfile(char(savepath),'02_preprocessed/','amica');
                mod = loadmodout12(icapath);
                
                %apply ICA weights to data
                EEGcontinuous.icasphere = mod.S;
                EEGcontinuous.icaweights = mod.W;
                EEGcontinuous = eeg_checkset(EEGcontinuous);
                
                %find corresponding sensor indices
                %for k=1:numel(EEGcontinuous.chanlocs)
                %   str{k}=EEGcontinuous.chanlocs(k).labels;
                %end
                
                
                eeglab redraw
                
                %find corresponding sensor indices
                for k=1:numel(EEG.chanlocs), str{k}=EEG.chanlocs(k).labels; end
                
                EEG.nbchan = size(EEG.data,1);
                
                
                %%% MARK COMPONENTS TO REJECT %%%
                if autorun == 0 % manual mode
                    pop_selectcomps(EEG, [1:length(EEG.icachansind)] );
                    
                    choice=menu('Components marked?','yes','no - interrupt');
                    if choice==2
                        error('Interrupt requested')
                    end
                    
                    comps_to_rej = find(EEG.reject.gcompreject);
                    
                    
                else
                    tmp = load(fullfile(char(savepath),'02_preprocessed/'), sprintf('5_%s_ICAcleancont.mat',setname),'comps_to_rej');
                    comps_to_rej = tmp.comps_to_rej;
                    
                end
                
                %remove marked components from CONTINOUS data set
                %we keep the EEG continous so we can later on epoch (or not) in any way we want
                EEG = pop_subcomp( EEGcontinuous, comps_to_rej, 0);
                
                %save file
                EEG.preprocessing = [EEG.preprocessing 'ICACleanedcont,'];
                EEG = pop_editset(EEG, 'setname', sprintf('5_%s_ICAcleancont',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('5_%s_ICAcleancont',setname),'filepath',fullfile(char(savepath),'02_preprocessed/'));
                save(fullfile(char(savepath),'02_preprocessed/',sprintf('5_%s_ICAcleancont.mat',setname)),'comps_to_rej');

                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                %%
            case 7 %interpolate, reference
                %load data
                EEG = pop_loadset(sprintf('5_%s_ICAcleancont.set',setname),fullfile(char(savepath),'02_preprocessed/'));
                
                %%% INTERPOLATE CHANNELS %%%
                
                %load original chanlocs and target channels
                load(fullfile(char(savepath),'02_preprocessed/',sprintf('2_%s_channelrejTriggersXensor.mat',setname)));
                
                %which channels should NOT be interpolated?
                alldel = {'CzAmp2' 'BIP1' 'BIP2' 'BIP3' 'BIP4' 'BIP5' 'BIP6' 'BIP7' 'BIP8' 'AUX1' 'AUX2' 'AUX3' 'AUX4' 'AUX5' 'AUX6' 'AUX7' 'AUX8','TIME' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA' 'R_GAZE_X' 'R_GAZE_Y' 'R_AREA' 'INPUT' 'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'R-GAZE-X' 'R-GAZE-Y' 'R_AREA'};
                
                %to be sure to have same channel structure for all subjects 
                %we need to delete the VEOG (that one is not interpolated 
                %and might be deleted for some subjects but not for others) 
                %also to get the correct channellocations for topoplots 
                %later on (as empty channels will be ignored by function
                %'topoplot' which shifts all channeles by one
                EEG = pop_select(EEG, 'nochannel', {'VEOG'});
                
                
                % find which channels should be interpolated
                idxs = ~ismember({complete_chanlocs.labels},alldel);
                %interpolate missing channels
                EEG= pop_interp(EEG,complete_chanlocs(idxs),'spherical');
                
                EEG.preprocessing = [EEG.preprocessing ' channelInterpol'];
                
                %%% REREFERENCE %%% this step is ommited because we chose
                %%% Cz single reference earlier during processing
               % EEG = pop_reref( EEG, []); %Participants’ averages were then re-referenced to a common average reference. (Rossion & Caharel, 2011)
               % EEG.preprocessing = [EEG.preprocessing 'Rereference,'];
                
                %save files
                EEG = pop_editset(EEG, 'setname', sprintf('6_%s_Interp',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('6_%s_Interp',setname),'filepath',fullfile(char(savepath),'02_preprocessed/'));
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
        end
        
        test_preproc_finished;
        set(h,'Visible','on')
        if autorun == 0 % manual mode
            uiwait(h)
            selection = getappdata(h, 'selection');
        end
        
    end
    
    %% Catch an error if one occured
catch ME
    close(h)
    rethrow(ME)
end







