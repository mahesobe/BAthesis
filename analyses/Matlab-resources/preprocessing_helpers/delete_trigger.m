function [behav,EEG,triggernames] = delete_trigger(behav,EEG)

todelete=[];

%reshape values in delphi struct
behav.delphi.n170_cat=reshape(behav.delphi.n170_cat',size(behav.delphi.n170_cat,1)*size(behav.delphi.n170_cat,2),1);
behav.delphi.n170_picID=reshape(behav.delphi.n170_picID',size(behav.delphi.n170_picID,1)*size(behav.delphi.n170_picID,2),1);
            
%yield expected triggers from behavioral file and run rough compare
expected_triggers = reshape(behav.cat_ids',size(behav.cat_ids,1)*size(behav.cat_ids,2),1);
  
triggernames = unique(expected_triggers);
fprintf('TRIGGERS: %s\n',mat2str(triggernames))
        
x=size(EEG.event);
for t=1:max(x)
	EEG.event(t).type = deblank(EEG.event(t).type);
end
        
starts=find(strcmp({EEG.event.type},'200'));
if length(starts)>1
	delete_idx=max(starts);
	todelete=1:delete_idx-1;
	EEG.event(todelete)=[];
	EEG.deleted_events=todelete;
end
        
        
        
        
        
        
