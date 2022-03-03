data = guihandles(h);

set(data.step1text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step2text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step3text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step4text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step5text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step6text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step7text, 'String', '','enable','on','ForegroundColor',[0 0 0]);
set(data.step1, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step2, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step3, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step4, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step5, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step6, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);
set(data.step7, 'String', 'run','enable','on','ForegroundColor',[0 0 0]);


if exist(fullfile(filepath,'preprocessed',sprintf('2_%s_bandpass_resample_deblank.set',setname)),'file')
    set(data.step1text, 'String', 'already run');
    set(data.step1text, 'ForegroundColor',[1 0 0])
    set(data.step1,'string','rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed',sprintf('3_%s_channelrejTriggersXensor.set',setname)),'file')
    set(data.step2text, 'String', 'already run');
    set(data.step2text, 'ForegroundColor',[1 0 0])
    set(data.step2,'string','rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed',sprintf('4_%s_Clean.set',setname)),'file')
    set(data.step3text, 'String', 'already run');
    set(data.step3text, 'ForegroundColor',[1 0 0])
    set(data.step3,'string','rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed','amica','W'),'file')
    set(data.step4text, 'String', 'already run');
    set(data.step4text, 'ForegroundColor',[1 0 0])
    set(data.step4,'string','rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed',sprintf('5_%s_ICAEpoched.set',setname)),'file')
    set(data.step5text, 'String', 'already run');
    set(data.step5text, 'ForegroundColor',[1 0 0])
    set(data.step5, 'String', 'rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed',sprintf('6_%s_ICAcleancont.set',setname)),'file')
    set(data.step6text, 'String', 'already run');
    set(data.step6text, 'ForegroundColor',[1 0 0])
    set(data.step6, 'String', 'rerun','ForegroundColor','red','enable','on');
end

if exist(fullfile(filepath,'preprocessed',sprintf('7_%s_RerefInterp.set',setname)),'file')
    set(data.step7text, 'String', 'already run');
    set(data.step7text, 'ForegroundColor',[1 0 0])
    set(data.step7, 'String', 'rerun','ForegroundColor','red','enable','on');
end
