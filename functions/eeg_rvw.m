function eeg_rvw(iSubj)

global RUN;
switch RUN.dir.sess
    case 'enc', pre=RUN.enc.pre;  pdir='1';
    case 'ret', pre=RUN.ret.pre;  pdir='2';
    case 'err', pre=RUN.ret.pre;  pdir='3';
    otherwise, pre=RUN.pre; pdir='';
end

reject={'rejthreshE','rejconstE'}; RGB_mat=jet(2);

sdisp(RUN.dir.inc_subjects{iSubj},1);
%=================================%
% Check channel data
%=================================%
load_dir=fullfile(RUN.dir.pre,RUN.dir.subjects{iSubj},pdir,filesep);
post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{4}];
EEG=pop_loadset('filename',post_save,'filepath',load_dir);

load(fullfile(RUN.dir.pro,RUN.dir.inc_subjects{iSubj},...
    pdir,['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]),'data'); 
% Create win reject matrix
r=1;

Window=625;

Time_mat=[1:Window:Window*(EEG.trials+1);Window:Window:Window*(EEG.trials+1)]';
for ii=1:length(reject)
    for jj=1:length(EEG.reject.(reject{ii}(1:end-1)))
        if EEG.reject.(reject{ii}(1:end-1))(jj)==1
           % Found reject trial
           reject_mat(r,:)=[Time_mat(jj,:) RGB_mat(ii,:) ...
               EEG.reject.(reject{ii})(:,jj)'];
           rtype(r)=ii;
           r=r+1;
        end
    end
end

eegplot(EEG.data,'winrej',reject_mat((or(rtype==1,rtype==2)),:));


