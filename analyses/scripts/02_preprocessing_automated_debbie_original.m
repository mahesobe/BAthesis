% This script allows you to preprocess your EEG data
% semiautomatically.
%
% It contains the following steps:
% 1 - loading, filtering, downsampling, deblanking
% 2 - behavioral data, channel position, channel cleaning
% 3 - cleaning continous data
% 4 - calculate ICA
% 5 - apply ICA weights and temporary epoching
% 6 - show ICA components, clean them (from continous data)
% 7 - interpolation of missing channels, rereferencing to average
%
% Every step saves the resulting EEG file on store.
% In the end you will have a clean, continous data set with all channels
% and average refrerence.
%
%
% The script developed over time, but mostly it was written by @tkietzma (automation, GUI)
% @sitimm (adjusting, behavioral parts), @agert (adjusting, making it ready for EEG training)
%% to restart everything
% do not run when running it automatically
% clear variables
% close all;
% clc;
% restoredefaultpath
%% Load necessary code
cd /net/store/nbp/projects/wd_ride_village/processedData/village %folder to save matlab steps
savepath = "/net/store/nbp/projects/wd_ride_village/processedData/village";

% add all necesarry tools
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/preprocessing_helpers');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/NoiseTools');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2021.0');
remove_from_path('eeglab');  % remove old eeglab
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0/plugins/xdfimport1.18/xdf');

eeglab;

% Where is your data stored?
basepath='/net/store/nbp/projects/wd_ride_village/recordedData/wd_village';


%%% LOAD RAW DATA %%%
%select data file to work on

%cntpath = fullfile(basepath,sprintf('VP%u',sub),sprintf('EEGTraining_subj%u.cnt',sub));
% make sure to change that before proper analysis --> file names need to be
% changed
cntpath = fullfile(basepath,'04_v_020621.xdf');
if exist(cntpath,'file')
    [filepath,filename,ext] = fileparts(cntpath); % define variables
    % to keep the same as the manual selection procedure
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
    sub = 4;
    
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
                mkdir(fullfile(savepath,'preprocessed/'))
%                 command = ['chgrp -hR nbp ' fullfile(savepath,'preprocessed/')];
%                 system(command);
%                 command = ['chmod 771 -R ' fullfile(filepath,'preprocessed/')];
%                 system(command);
                  % not sure yet why this part is not working --> check
                  % again
%                 
                
                % load data into Mobilab 
                allDataStreams = dataSourceXDF(filename,'MobiLabData');
                EEG = allDataStreams.item{9};
                EEG = allDataStreams.export2eeglab( 9 )
                
                path = fullfile(savepath,'preprocessed/');
                EEG = pop_saveset( EEG, 'filename', sprintf('0_%s_raw',setname), 'filepath', char(fullfile(savepath,'preprocessed/')));

                
                EEG = pop_editset(EEG, 'setname', sprintf('0_%s_raw',setname));
                EEG = eeg_checkset(EEG)

                
                EEG.preprocessing = 'Raw,'; %every step gets stored in .preprocessing so you can see what happened to the data set
                
%                 % load the before mentioned .set into eeglab
%                 sub = 'pilot1';
%                 filename = 'new_pilot1.set';
%                 filepath = 'E:/Uni/Bachelor';
%                 EEG1 = pop_loadset(filename, filepath);
                %
%                 EEG.setname = filename;
%                 setname = EEG.setname;
%                 EEG = eeg_checkset(EEG);
                
                % make EEG file that is untouched for later use
                trEEG=EEG;
                
                
                %%% HighPass FILTER %%%
                % filter out low frequencies
                if strcmp(filtType,'acausal')
                    EEG = pop_eegfiltnew(EEG, filtFreq, []);
                elseif strcmp(filtType,'causal')
                    EEG = pop_eegfiltnew(EEG, filtFreq,[],[],[],[],0,1);
                end
                
                EEG.preprocessing = [EEG.preprocessing 'Highpass' num2str(filtFreq) ', '];
                
                
                %%% DOWNSAMPLING %%%
                % lower sampling rate from default (1024Hz) to 512Hz
                EEG = pop_resample(EEG, 512);
                EEG.preprocessing = [EEG.preprocessing 'Resampled,'];
                trEEG = pop_resample(trEEG, 512);
                
                
                % Save the data
                EEG = pop_editset(EEG, 'setname', sprintf('2_%s_bandpass_resample_deblank',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('2_%s_bandpass_resample_deblank',setname),'filepath',char(fullfile(savepath,'preprocessed/')));
                
                trEEG = pop_saveset(trEEG, 'filename',sprintf('%s_triggers',setname),'filepath',char(fullfile(savepath,'preprocessed/')));
                
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 2 %behav, position, channelselect
                
                fprintf('Loading filtered and resampled files...\n')
                
                %load files
                EEG = pop_loadset(sprintf('2_%s_bandpass_resample_deblank.set',setname),fullfile(filepath,'preprocessed/' ));
                trEEG = pop_loadset(sprintf('%s_triggers.set',setname),fullfile(filepath,'preprocessed/' ));
                
                

                
                %% Channel positions
                % if you used the Xensor
                if exist(fullfile(basepath,sprintf('VP%u',sub),'*.elc'), 'file') == 2
                    [xfilename, xfilepath]=uigetfile(fullfile(basepath,sprintf('VP%u',sub),'*.elc'));
                    EEG = load_xensor(EEG,fullfile(xfilepath,xfilename));
                else %use default positions
                    EEG=pop_chanedit(EEG, 'lookup','/net/store/nbp/projects/EEG_Training/Analysis/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
                end
                EEG.preprocessing = [EEG.preprocessing 'ChanPos,'];
                EEG_tmp=EEG;
                
                
                %%% CHANNEL CLEANING %%%
                % which channels do not contain EEG data by default?
                alldel = {'BIP2' 'BIP3' 'BIP4' 'BIP5' 'BIP6' 'BIP7' 'BIP8' 'AUX1' 'AUX2' 'AUX3' 'AUX4' 'AUX5' 'AUX6' 'AUX7' 'CZ' 'IZ' 'TIME' 'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'R-GAZE-X' 'R-GAZE-Y' 'R-AREA' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA', 'R_GAZE_X' 'R_GAZE_Y' 'R_AREA' 'INPUT'};
                
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
                if exist(fullfile(filepath,'preprocessed/',sprintf('3_%s_channelrejTriggersXensor.mat',setname)))
                    tmp=load(fullfile(filepath,'preprocessed',sprintf('3_%s_channelrejTriggersXensor.mat',setname)));
                    
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
                
                %save file
                EEG = pop_editset(EEG, 'setname', sprintf('3_%s_channelrejTriggersXensor',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('3_%s_channelrejTriggersXensor',setname),'filepath',fullfile(filepath,'preprocessed/' ));
                
                save(fullfile(filepath,'preprocessed/',sprintf('3_%s_channelrejTriggersXensor.mat',setname)),'chan_del','trigger_names','targetchannels','complete_chanlocs');
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
                
                
            case 3 %data cleaning
                fprintf('Loading channel cleaned files...\n')
                
                
                %load files
                EEG = pop_loadset(sprintf('3_%s_channelrejTriggersXensor.set',setname),fullfile(filepath,'preprocessed/'));
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
                if exist(fullfile(filepath,'preprocessed',sprintf('4_%s_cleaningTimes.mat',setname)),'file')
                    load(fullfile(filepath,'preprocessed',sprintf('4_%s_cleaningTimes.mat',setname)),'tmprej','rej');
                    
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
                EEG = pop_editset(EEG, 'setname', sprintf('4_%s_Clean',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('4_%s_Clean',setname),'filepath',fullfile(filepath,'preprocessed/'));
                
                save(fullfile(filepath,'preprocessed/',sprintf('4_%s_cleaningTimes.mat',setname)),'tmprej','rej');
                
                
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
                
                if exist(fullfile(filepath,'preprocessed/','amica'))
                    fprintf('ICA folder found')
                    %wait for to confirm that ICA shall be rerun
                    choice = menu('Do you really want to rerun the ICA?','yes');
                end
                
                
                %load files
                EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed/'));
                
                
                %to exclude luminance sensor channels from ICA (data not lost, but excluded from analysis)
                auxdel = find(cellfun(@(x)~isempty(x),strfind({EEG.chanlocs(:).labels},'AUX')));
                EEG = pop_select(EEG, 'nochannel', auxdel);
                
                
                %for the ICA computation, highpass-filter data at 2Hz to get better
                %Eye-Components
                EEG2Hz = pop_eegfiltnew(EEG, 2, []);   % highpass
                
                
                %%% ICA %%%
                addpath('./amica')
                mkdir(fullfile(filepath,'preprocessed/','amica'))
                permission_cleanup(filepath);
                outDir = fullfile(filepath,'preprocessed/','amica');
                
                
                %run ICA on grid
                runamica12(EEG2Hz.data,'outdir',outDir,'use_queue','all.q','qsubname',['ica_VP' num2str(sub)]);
                
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 5 %Apply ICA weights and Epoch
                fprintf('Loading data cleaned files and ICA information...\n')
                
                addpath([fullfile(filepath,'preprocessed/') '/amica'])
                
                
                %load files
                EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed/' ));
                load(fullfile(filepath,'preprocessed/',sprintf('3_%s_channelrejTriggersXensor.mat',setname)))
                
                
                %load ICA results
                icapath = fullfile(filepath,'preprocessed/','amica');
                mod = loadmodout12(icapath);
                
                
                %%% APPLY ICA WEIGHTS TO DATA %%%
                EEG.icasphere = mod.S;
                EEG.icaweights = mod.W;
                EEG = eeg_checkset(EEG);
                EEG.preprocessing = [EEG.preprocessing 'AMICA,'];
                
                
                %%% EPOCH DATA %%%
                % this epoching is done only for the vizualization later on
                % we want to keep the actual data set continous so that we are not restricted by this step later on
                epoch_trigg={'1','2','3','4','5'};
                window     =[-0.1 0.3];
                baseline   =[-100 0];
                
                [EEG,indices] = pop_epoch( EEG, epoch_trigg, window, 'epochinfo', 'yes');
                EEG.orig_indices = indices;
                
                EEG = pop_rmbase( EEG, baseline);
                EEG.preprocessing = [EEG.preprocessing 'Epoched,'];
                
                
                %save files
                EEG = pop_editset(EEG, 'setname', sprintf('5_%s_ICAEpoched',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('5_%s_ICAEpoched',setname),'filepath',fullfile(filepath,'preprocessed/' ));
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                
            case 6 %component cleaning
                comps_to_rej = [];
                
                %load files
                EEG = pop_loadset(sprintf('5_%s_ICAEpoched.set',setname),fullfile(filepath,'preprocessed/'));
                EEGcontinuous= pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed/' ));
                
                icapath = fullfile(filepath,'preprocessed/','amica');
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
                    tmp = load(fullfile(filepath,'preprocessed/',sprintf('6_%s_ICAcleancont.mat',setname)),'comps_to_rej');
                    comps_to_rej = tmp.comps_to_rej;
                    
                end
                
                %remove marked components from CONTINOUS data set
                %we keep the EEG continous so we can later on epoch (or not) in any way we want
                EEG = pop_subcomp( EEGcontinuous, comps_to_rej, 0);
                
                %save file
                EEG.preprocessing = [EEG.preprocessing 'ICACleanedcont,'];
                EEG = pop_editset(EEG, 'setname', sprintf('6_%s_ICAcleancont',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('6_%s_ICAcleancont',setname),'filepath',fullfile(filepath,'preprocessed/'));
                save(fullfile(filepath,'preprocessed/',sprintf('6_%s_ICAcleancont.mat',setname)),'comps_to_rej');
                
                fprintf('permission cleanup\n')
                permission_cleanup(filepath);
                fprintf('done\n')
                test_preproc_finished;
                fprintf('gui update done\n')
                
                %%
            case 7 %interpolate, reference
                %load data
                EEG = pop_loadset(sprintf('6_%s_ICAcleancont.set',setname),fullfile(filepath,'preprocessed/'));
                
                %%% INTERPOLATE CHANNELS %%%
                
                %load original chanlocs and target channels
                load(fullfile(filepath,'preprocessed/',sprintf('3_%s_channelrejTriggersXensor.mat',setname)));
                
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
                
                %%% REREFERENCE %%%
                EEG = pop_reref( EEG, []); %Participants’ averages were then re-referenced to a common average reference. (Rossion & Caharel, 2011)
                EEG.preprocessing = [EEG.preprocessing 'Rereference,'];
                
                %save files
                EEG = pop_editset(EEG, 'setname', sprintf('7_%s_RerefInterp',setname));
                EEG = pop_saveset(EEG, 'filename',sprintf('7_%s_RerefInterp',setname),'filepath',fullfile(filepath,'preprocessed/'));
                
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







