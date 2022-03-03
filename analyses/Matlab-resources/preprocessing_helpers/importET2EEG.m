% synchronize EyeLink and EEG
% uses EYEEEG toolbox by Diemiegen [1]

% [1] Dimigen, O., Sommer, W., Hohlfeld, A., Jacobs, A., & Kliegl, R. (2011). 
% Coregistration of eye movements and EEG in natural reading: Analyses & Review. 
% Journal of Experimenta Psychology: General, 140 (4), 552-572 
% http://www2.hu-berlin.de/eyetracking-eeg

% @sitimm, @agert

%This is an example script from IntoTheWild. You have to adjust the pop_importeyetracker line
% with your own start, stop and triggers of interest.


function EEG = importET2EEG(EEG,path,sub)

% functions need EEG.urevent - fill with data from EEG.event, as urevent
% ist still empty (raw data)EEG = pop_importeyetracker(EEG, [path '/asc2mat.mat'] ,[200 255]);
EEG.urevent=EEG.event;

for e=1:length(EEG.event)
    EEG.event(e).duration=0;
    EEG.event(e).urevent=e;
end


% here we need to match a regular expression in the code
% message of sending sth via parallel port -> event data

    ET = parseeyelink([path '/' 'subj' num2str(sub) '.asc'],[path '/asc2mat.mat'],'!CMD . write_ioport 0x378');
    
% actually importing the data into the EEG

    EEG = pop_importeyetracker(EEG, [path '/asc2mat.mat'] ,[200 255],[1:5] ,{'TIME' 'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'INPUT'},1,1,0,0,20);

%to plot and save the figure of the synchronozation output, change last 0
%in pop_importeyetracker into a 1

