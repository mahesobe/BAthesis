% This script adds trial information like the number of the stimulus picture or the subject number into
% the EEG structure.
% Here you would also add things like correctnes of response, reaction time, etc,
% Furthermore, it compares the triggers that should have been sent with those that were sent, and deletes all 
% of those that were not correct (=ghose trigggers). If more than 7 are found, most likely you did sth wrong.

function [behav,EEG,trEEG,triggernames] = insert_behav_WLFO (behav,EEG,trEEG)

todelete=[];

%reshape values in delphi struct
behav.delphi.n170_cat=reshape(behav.delphi.n170_cat',size(behav.delphi.n170_cat,1)*size(behav.delphi.n170_cat,2),1);
behav.delphi.n170_picID=reshape(behav.delphi.n170_picID',size(behav.delphi.n170_picID,1)*size(behav.delphi.n170_picID,2),1);

%yield expected triggers from behavioral file and run rough compare
expected_triggers = reshape(behav.cat_ids',size(behav.cat_ids,1)*size(behav.cat_ids,2),1);
triggernames = unique(expected_triggers);


fprintf('TRIGGERS: %s\n',mat2str(triggernames))

trial_count=1;

for eventid=2:length(EEG.event)
    if length(EEG.event(eventid).type)<=5
        
        if isempty(find(triggernames==str2num(EEG.event(eventid).type)))==0
            if str2num(EEG.event(eventid).type)~=expected_triggers(trial_count)
                fprintf('Deleting %u:\t%s\n',eventid,EEG.event(eventid).type);
                todelete = [todelete eventid];
            else
                
                EEG.event(eventid).id          = behav.delphi.subjectID(trial_count);
                EEG.event(eventid).stimCat     = behav.delphi.n170_cat(trial_count);
                EEG.event(eventid).picID       = behav.delphi.n170_picID(trial_count);
                EEG.event(eventid).origLatency = EEG.event(eventid).latency;
                EEG.event(eventid).trialnum    = trial_count;
                
                trial_count=trial_count+1;
                
            end
        end
    end
end


if length(todelete) > 7
    error('Error! current subject has too many trials to be deleted!')
end
EEG.deleted_events=todelete;
EEG.event(todelete)=[];

EEG.preprocessing = [EEG.preprocessing,'Behavioral,'];

%change events here



for b = 1 : length(trEEG)
    trEEG.event(b).id          = 0;
    trEEG.event(b).stimCat     = 0;
    trEEG.event(b).picID       = 0;
    trEEG.event(b).origLatency = 0;
end

triggernames=num2cell(triggernames);
