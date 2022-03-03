% CAUTION! Example script from IntoTheWild
% This function takes Lab data and calculates EEG.event.humanface etc from
% the EEG triggers.
function EEG=addbehavtoLab(EEG)

%% find current and prev info
HF=sort([find(strcmp({EEG.event.type},'1')),find(strcmp({EEG.event.type},'2')),find(strcmp({EEG.event.type},'3'))]);
NH=sort([find(strcmp({EEG.event.type},'4')),find(strcmp({EEG.event.type},'5'))]);

evts(1,:)=sort([HF NH]);

currcat={EEG.event(evts).type};
evts(2,:)=cellfun(@str2num,currcat);

prevcat=[6 evts(2,1:end-1)];
for b=1:7
    prevcat(160*b+1)=6;
end
evts(3,:)=prevcat;

%% initialize all fields
numevts=length(EEG.event);

humanface=zeros(1,numevts);
prevhumanface=zeros(1,numevts);
prevprevhumanface=zeros(1,numevts);
humanhead=zeros(1,numevts);
prevhumanhead=zeros(1,numevts);
prevprevhumanhead=zeros(1,numevts);
nonhuman=zeros(1,numevts);
prevnonhuman=zeros(1,numevts);
prevprevnonhuman=zeros(1,numevts);
none=zeros(1,numevts);
prevnone=zeros(1,numevts);
prevprevnone=zeros(1,numevts);
overlapping=zeros(1,numevts);
prevoverlapping=zeros(1,numevts);
prevprevoverlapping=zeros(1,numevts);
outside=zeros(1,numevts);
prevoutside=zeros(1,numevts);
prevprevoutside=zeros(1,numevts);
samebox=zeros(1,numevts);
prevsamebox=zeros(1,numevts);
sameboxn2=zeros(1,numevts);
notwithintrial=zeros(1,numevts);

%% add into vectors
humanface(HF)=1;
none(NH)=1;
prevhumanface(evts(1,evts(3,:)<4))=1;
a=evts(3,:)>3;
b=evts(3,:)<6;
ab=a+b;
prevnone(evts(1,ab==2))=1;

%% add into EEG structure
for ne=1:numevts
    EEG.event(ne).humanface=humanface(ne);
    EEG.event(ne).prevhumanface=prevhumanface(ne);
    EEG.event(ne).prevprevhumanface=prevprevhumanface(ne);
    EEG.event(ne).humanhead=humanhead(ne);
    EEG.event(ne).prevhumanhead=prevhumanhead(ne);
    EEG.event(ne).nonhuman=nonhuman(ne);
    EEG.event(ne).prevnonhuman=prevnonhuman(ne);
    EEG.event(ne).prevprevnonhuman=prevprevnonhuman(ne);
    EEG.event(ne).none=none(ne);
    EEG.event(ne).prevnone=prevnone(ne);
    EEG.event(ne).prevprevnone=prevprevnone(ne);
    EEG.event(ne).overlapping=overlapping(ne);
    EEG.event(ne).prevoverlapping=prevoverlapping(ne);
    EEG.event(ne).prevprevoverlapping=prevprevoverlapping(ne);
    EEG.event(ne).outside=outside(ne);
    EEG.event(ne).prevoutside=prevoutside(ne);
    EEG.event(ne).prevprevoutside=prevprevoutside(ne);
    EEG.event(ne).samebox=samebox(ne);
    EEG.event(ne).prevsamebox=prevsamebox(ne);
    EEG.event(ne).sameboxn2=sameboxn2(ne);
    EEG.event(ne).notwithintrial=notwithintrial(ne);
    
end

