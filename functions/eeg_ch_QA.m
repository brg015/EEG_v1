function EEG=eeg_ch_QA(EEG,pre,iSubj,sub_dir,skip)

global RUN; close all;
%=========================================================================%
% Collect Channel Measures
%=========================================================================%
pre.thresh.epc=RUN.set.thresh_epc;

ch_set=1:EEG.nbchan;
if skip==0

    % Reject trials from all electrodes but they eye
    [EEG,~] = pop_eegthresh(EEG,1,ch_set,pre.thresh_dta(1),pre.thresh_dta(2),...
        pre.thresh.epc(1),pre.thresh.epc(2),1,0);

    local_thresh=3;  % 3 sigma events
    global_thresh=3; % 3 sigma events
    [EEG, ~, ~, ~]=pop_jointprob(EEG,1,ch_set,local_thresh,global_thresh,1,0,0);

    % Looking at 1 second
    max_slope=50;
    minR=0.3;
    EEG = pop_rejtrend(EEG,1,ch_set,RUN.set.sf,max_slope,minR,1,0,0);

    % Muscle
    dB_thresh=[-100 25];
    fq_thresh=[20 40];
    [EEG,spec_indm] = pop_rejspec(EEG,1,'threshold',dB_thresh,...
        'freqlimits',fq_thresh,'eegplotreject',0,'eegplotplotallrej',1);

    % Eye
    dB_thresh=[-50 50];
    fq_thresh=[0 2];
    [EEG_eye,spec_inde] = pop_rejspec(EEG,1,'threshold',dB_thresh,...
        'freqlimits',fq_thresh,'eegplotreject',0,'eegplotplotallrej',1);

     % Default
    dB_thresh=[-25 25];
    fq_thresh=[0 50];ls
    [EEG_def,~] = pop_rejspec(EEG,1,'threshold',dB_thresh,...
        'freqlimits',fq_thresh,'eegplotreject',0,'eegplotplotallrej',1);

    % Combine eye and muscle data
    EEG.reject.rej_eye_freq=EEG_eye.reject.rejfreq;
    EEG.reject.rej_eye_freqE=EEG_eye.reject.rejfreqE;
    EEG.reject.rej_def_freq=EEG_def.reject.rejfreq;
    EEG.reject.rej_def_freqE=EEG_def.reject.rejfreqE;
    clear EEG_def EEG_eye;
end
%=========================================================================%
% Save that sexy output
%=========================================================================%
save_dir=fullfile(RUN.dir.QAL,sub_dir);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

reject={'rejjpE','rejthreshE','rejconstE' 'rejfreqE','rej_eye_freqE' 'rej_def_freqE'};
% Remove eye electrode details.
if skip==1
     reject={'rejthreshE' 'rejconstE'};
    % Reject trials from all electrodes but they eye
    [EEG,~] = pop_eegthresh(EEG,1,setdiff(ch_set,31),pre.thresh_dta(1),pre.thresh_dta(2),...
        pre.thresh.epc(1),pre.thresh.epc(2),1,0);
    
    max_slope=50;
    minR=0.3;
    EEG = pop_rejtrend(EEG,1,setdiff(ch_set,31),RUN.set.sf,max_slope,minR,1,0,0);
end

r=1; Window=625; RGB_mat=jet(6); 
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
%===================%
% Rejections
%===================%
for ii=1:length(reject)
    indx=(rtype==ii);
    reject_count(ii,:)=sum(reject_mat(indx,6:end),1);
    if skip==0
        reject_z(ii,:)=zscore(reject_count(ii,:));
    else
        reject_z(ii,:)=reject_count(ii,:)/size(reject_mat,1);
    end
    figure(1); plot(reject_count(ii,:),'color',RGB_mat(ii,:),'linewidth',3); hold on;
    figure(2); plot(reject_z(ii,:),'color',RGB_mat(ii,:),'linewidth',3); hold on;  
end

if skip==0
    figure(1); set(gcf,'position',[0 0 1280 1024]);
    plot(mean(reject_count),'k','linewidth',3); legend(reject); grid;
    set(gca,'Color',[0.6 0.6 0.6]); 
    export_fig(fullfile(save_dir,[RUN.dir.subjects{iSubj} '_rej_count.png']));
    close(1);

    figure(2); set(gcf,'position',[0 0 1280 1024]);
    plot(mean(reject_z),'k','linewidth',3); legend(reject); grid;
    set(gca,'Color',[0.6 0.6 0.6]);
    export_fig(fullfile(save_dir,[RUN.dir.subjects{iSubj} '_rej_z.png']));
    close(2);

    figure(3); set(gcf,'position',[0 0 1280 1024]);
    plot(mean(reject_z,1),'linewidth',3); grid;
    set(gca,'Color',[0.6 0.6 0.6]);
    export_fig(fullfile(save_dir,[RUN.dir.subjects{iSubj} '_mean_z.png']));
    close(3);
else
    figure(1); set(gcf,'position',[0 0 1280 1024]);
    plot(mean(reject_count),'k','linewidth',3); legend(reject); grid;
    set(gca,'Color',[0.6 0.6 0.6]); 
    export_fig(fullfile(save_dir,[RUN.dir.subjects{iSubj} '_rej_count.png']));
    close(1);
    
    figure(2); set(gcf,'position',[0 0 1280 1024]);
    plot(mean(reject_z),'k','linewidth',3); legend(reject); grid;
    set(gca,'Color',[0.6 0.6 0.6]);
    export_fig(fullfile(save_dir,[RUN.dir.subjects{iSubj} '_rej_per.png']));
    close(2);
end

%================%
% Heat maps
%================%
rej_crit=reject; g=[];
for ii=1:length(rej_crit),
    rej_save=fullfile(save_dir,[RUN.dir.subjects{iSubj} '_' rej_crit{ii} '.png']);
    eeg_pdm(EEG.reject.(rej_crit{ii}),'','',rej_save);
    if ii>1,
        g=g+EEG.reject.(rej_crit{ii});
    else
        g=EEG.reject.(rej_crit{ii});
    end
end
rej_save=fullfile(save_dir,[RUN.dir.subjects{iSubj} '_all.png']);
eeg_pdm(g,'','',rej_save);
    
SL.subj_per(iSubj)=sum(or(rtype==2,rtype==3))/length(rtype);



