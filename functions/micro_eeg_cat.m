function data=micro_eeg_cat(subj_str)
global RUN;
skip_stuff=0;

pro_dir=fullfile(RUN.dir.pro,subj_str);

if exist(fullfile(pro_dir,'original'),'dir'), skip_stuff=1; end

% if folder already exists, don't need to make stuff and move things into
% it.
if skip_stuff==0
    mkdir(fullfile(pro_dir,'original'));

    movefile(fullfile(pro_dir,'2'),fullfile(pro_dir,'original','2'));
    movefile(fullfile(pro_dir,'3'),fullfile(pro_dir,'original','3'));
end

cd(fullfile(pro_dir,'original'));
wrk_dir=cd;
[save_dir,~,~]=fileparts(wrk_dir);

I=findstr('ret',RUN.save_str{5});
f2_file=RUN.save_str{5};
f2_file(I:I+2)='err';
f1=fullfile(wrk_dir,'2',['subj_' subj_str '_' RUN.save_str{5}]);
f2=fullfile(wrk_dir,'3',['subj_' subj_str '_' f2_file]);
save_file=fullfile(pro_dir,'2',RUN.save_str{5});

load(f1); d1=data; clear data;
load(f2); d2=data; clear data;

% Use d1 as template 
data.label=d1.label;
data.fsample=d1.fsample;
data.elec=d1.elec;
data.cfg=d1.cfg;
% Ignore EEGepoch & EEGurevent as they're not used by ft
% Trial Struct is however needed for correct coding
data.trial=d1.trial;
data.time=d1.time;
data.trialStruct=d1.trialStruct;

N=length(data.time);
for ii=1:length(d2.time)
    data.trial{N+ii}=d2.trial{ii};
    data.time{N+ii}=d2.time{ii};
    data.trialStruct.trial(N+ii)=d2.trialStruct.trial(ii);
    data.trialStruct.phase(N+ii)=d2.trialStruct.phase(ii);
    data.trialStruct.event(N+ii)=d2.trialStruct.event(ii);
    data.trialStruct.rProfile(N+ii,:)=d2.trialStruct.rProfile(ii,:);
    data.trialStruct.rejthresh(N+ii)=d2.trialStruct.rejthresh(ii);
end
data.trialStruct.eProfile=[d1.trialStruct.eProfile;d2.trialStruct.eProfile];
data.trialStruct.descrip=d1.trialStruct.descrip;

if ~exist(fullfile(save_dir,'2'),'dir'),mkdir_tree(fullfile(save_dir,'2')); end
save(save_file,'data');
