function eeg_adv()

global RUN;
switch RUN.dir.sess
    case 'enc', pdir='1'; pre.baseline=RUN.enc.pre.baseline;
    case 'ret', pdir='2'; pre.baseline=RUN.ret.pre.baseline;
    case 'err', pdir='3'; pre.baseline=RUN.ret.pre.baseline;
    otherwise, pdir='';
end

[~,u_save_str,~]=fileparts(RUN.save_str{5});
% Only need this as a template...
gcohs_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(sum(RUN.dir.analyze)) '_cohs.mat']);

%=========================================================================%    
%% RUN COVARATION
%=========================================================================%
if RUN.adv.cov.on==1
%=========================================================================%    
% Define model inputs       
covary_model='COCA';
model_name='COCA_v8';
RUN.adv.model_name=model_name;
F_method='v2';
% Freq Analysis
switch F_method
    case ''
        cfg_freq.method='mtmconvol';
        cfg_freq.output='powandcsd';
        cfg_freq.taper='hanning';
        cfg_freq.channel='all';
        cfg_freq.keeptrials='yes';
        cfg_freq.keeptapers='no';
        cfg_freq.foi=4:2:30;
        cfg_freq.t_ftimwin=6./cfg_freq.foi; % 7 cycles per window
        cfg_freq.toi=-1:0.1:2.5;
        for ii=1:length(cfg_freq.foi)
            Bands{ii}=[cfg_freq.foi(ii) cfg_freq.foi(ii)];
            Band_names{ii}=['Freq_' num2str(cfg_freq.foi(ii))];
        end
    case 'v2'
        cfg_freq.method='mtmconvol';
        cfg_freq.output='powandcsd';
        cfg_freq.taper='hanning';
        cfg_freq.channel='all';
        cfg_freq.keeptrials='yes';
        cfg_freq.keeptapers='no';
        cfg_freq.foi=4:1:30;
        cfg_freq.t_ftimwin=6./cfg_freq.foi; % 7 cycles per window
        cfg_freq.toi=-1:0.1:2.5;
        for ii=1:length(cfg_freq.foi)
            Bands{ii}=[cfg_freq.foi(ii) cfg_freq.foi(ii)];
            Band_names{ii}=['Freq_' num2str(cfg_freq.foi(ii))];
        end
end
%=========================================================================%    
smooth_data=1;
switch model_name
    case 'COCA_v2'
        glm_mod.model_coef={'Int','Freq' 'WordM' 'PictM' 'Maxn' 'Trial' 'LFreq' 'Freqaxn' 'Size'};
        for ii=1:11, glm_mod.model_coef{end+1}=['ctg_' num2str(ii)]; end
        glm_mod.model_ind=0;
        glm_mod.mem=1;
        glm_mod.ds=0;
    case 'COCA_v3'
        glm_mod.model_coef={'Int','Freq' 'WordM' 'PictM' 'Maxn' 'Trial' 'LFreq' 'Freqaxn'};
        glm_mod.model_ind=1;
        glm_mod.mem=1;
        glm_mod.ds=1;
        glm_mod.freq=0;
    case 'COCA_v4'
        glm_mod.model_coef={'Int','Freq' 'WordM' 'PictM' 'Maxn' 'Trial' 'LFreq' 'Freqaxn'};
        glm_mod.model_ind=0;
        glm_mod.mem=1;
        glm_mod.ds=1;
        glm_mod.freq=0;
    case 'COCA_v5'
        glm_mod.model_coef={'Int','WordM' 'PictM' 'Maxn'};
        glm_mod.model_ind=1;
        glm_mod.mem=0;
        glm_mod.ds=1;
        glm_mod.freq=0;
    case 'COCA_v6'
        glm_mod.model_coef={'Int','WordM' 'Freq' 'LFreq' 'Freqaxn'};
        glm_mod.model_ind=0;
        glm_mod.mem=1;
        glm_mod.ds=1;
        glm_mod.freq=1;
   case 'COCA_v7'
        glm_mod.model_coef={'Int','WordM' 'Freq' 'LFreq' 'Freqaxn'};
        glm_mod.model_ind=0;
        glm_mod.mem=0;
        glm_mod.ds=1;
        glm_mod.freq=1;
    case 'COCA_v8'
        % Actually has interaxn terms
        glm_mod.model_coef={'Int' 'WordM' 'Freq' 'LFreq', ...
            'FreqXWordM' 'LFreqXWordM' 'FreqXLFreq'};
        glm_mod.model_ind=0;
        glm_mod.mem=2;
        glm_mod.ds=0;
        glm_mod.freq=1;
        smooth_data=0;
    case 'COCA_v9'
        glm_mod.model_coef={'Int','WordM' 'Freq' 'LFreq'};
        glm_mod.model_ind=0;
        glm_mod.mem=2;
        glm_mod.ds=0;
        glm_mod.freq=1;
        smooth_data=0;
end  
        
switch covary_model
    case 'COCA'
        
        model_data=excel_reader(RUN.template.COCA);
        model_id=cell2num(model_data{1}.col);   % ID
        model_val=cell2num(model_data{7}.col);  % COCA
        model_dom=cell2num(model_data{13}.col);           % Domain
        model_cat=cell2num(model_data{14}.col);           % Category
        model_siz=cell2num(model_data{15}.col);           % Size
        model_lfq=cell2num(model_data{17}.col);           % Letter Freq
        model_ctc=cell2num(model_data{18}.col);           % Catch
        model_nlt=cell2num(model_data{21}.col);           % Letter Count
        % Some slight manipulations here
        model_val=log(model_val); % Distrib. more guassian
        % Then normalize the measures as well
        model_val=(model_val-mean(model_val))/std(model_val);
        model_lfq=zscore(model_lfq);   
end

% glm_mod
%   .model_ind (0/1)   => run models with independence
%   .mem (0/1)         => if 1, calc. subj. perf. from d'
%   .model_coef{}      => model coefficients names
%=========================================================================%
if RUN.adv.cov.calc==1
%=========================================================================%
% Load in dprime template
dprime_data=excel_reader(RUN.template.dprime);
% ID, HC W, HC P, LC W, LC P

sdisp('GLM',1);
for iSubj=1:length(RUN.dir.subjects)
    
    glm_save=fullfile(RUN.adv.save,[model_name '_' RUN.dir.subjects{iSubj} '.mat']);
    if (exist(glm_save,'file') && RUN.dir.overwrite==0), continue; end
%=========================================================================%
% 1) Load in subject data
%=========================================================================%
    sdisp(['Subject: ' RUN.dir.subjects{iSubj}],2);

    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,filesep);

    post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]; % Based on trials*
    load(fullfile(save_pro_dir,post_save),'data'); % var=data;
    
    LV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
    if (exist(LV,'file') && glm_mod.freq==1)
        load(LV); 
%=========================================================================%
% 1b) Preproc frequency data
%=========================================================================%
        % Bands
        % Band_names
        cfg_f.baselinetype='ztransformSession';
        cfg_f.baseline=[-1 2.5]; % Not used to z-transform
        cfg_f.mu=[];
        cfg_f.sigma=[];
        Fz=ft_freqbaseline(cfg_f,freqlock);

        % eeg_glm input var{trial}[Ch X Time];
        for ii=1:length(Band_names)
            for jj=1:size(Fz.powspctrm,1)
                Fget=and(Fz.freq>=Bands{ii}(1),Fz.freq<=Bands{ii}(2));
                % Average across wanted frequency band
                Fband{ii}{jj}=squeeze(mean(Fz.powspctrm(jj,:,Fget,:),3));
            end
        end
        % Fband{band}{trial}[Ch X Time]
    end
%=========================================================================%
% 2) Short micro to distribute trials
%=========================================================================%
    % Need to create a mem regressor & reject regressor as well, also want
    % to account for subject memory effects as well.
    %========================%
    % Create mem_val regressor
    %========================%
    mem_word_val=nan(1,size(data.trialStruct.rProfile,1));
    mem_pict_val=nan(1,size(data.trialStruct.rProfile,1));
    
    % Calculate weights & Normalize
    if (glm_mod.mem==1 || glm_mod.mem==2)
        subj_ind=strcmp(RUN.dir.subjects{iSubj},dprime_data{1}.col);
        W=[str2double(dprime_data{2}.col{subj_ind}),str2double(dprime_data{4}.col{subj_ind}) ...
            -str2double(dprime_data{4}.col{subj_ind}),-str2double(dprime_data{2}.col{subj_ind})];
        P=[str2double(dprime_data{3}.col{subj_ind}),str2double(dprime_data{5}.col{subj_ind}) ...
            -str2double(dprime_data{5}.col{subj_ind}),-str2double(dprime_data{3}.col{subj_ind})];
        if glm_mod.mem==1
            Pf=(P-mean(P))/std(P);
            Wf=(W-mean(W))/std(W); 
        else
            Wf=[W(1:2),-W(2),-W(2)];
            Pf=[P(1:2),-P(2),-P(2)];
        end
    else
        % These are defined wrong :( - not nay more:
        Wf=-[-1.5 -0.5 0.5 1.5];
        Pf=-Wf;
    end

    mem_word_val(and(data.trialStruct.rProfile(:,3),data.trialStruct.rProfile(:,4)))=Wf(1); 
    mem_word_val(and(~data.trialStruct.rProfile(:,3),data.trialStruct.rProfile(:,4)))=Wf(2);
    mem_word_val(and(~data.trialStruct.rProfile(:,3),data.trialStruct.rProfile(:,5)))=Wf(3);
    mem_word_val(and(data.trialStruct.rProfile(:,3),data.trialStruct.rProfile(:,5)))=Wf(4);
    mem_pict_val(and(data.trialStruct.rProfile(:,8),data.trialStruct.rProfile(:,9)))=Pf(1);
    mem_pict_val(and(~data.trialStruct.rProfile(:,8),data.trialStruct.rProfile(:,9)))=Pf(2);
    mem_pict_val(and(~data.trialStruct.rProfile(:,8),data.trialStruct.rProfile(:,10)))=Pf(3);
    mem_pict_val(and(data.trialStruct.rProfile(:,8),data.trialStruct.rProfile(:,10)))=Pf(4);
    % Create reject regressor
    reject_val=zeros(1,size(data.trialStruct.rProfile,1));
    reject_val(and(mem_word_val==0,mem_pict_val==0))=1;
    reject_val(logical(data.trialStruct.rejthresh))=1;
    % Update model_val regressor
    for ii=1:size(data.trialStruct.rProfile,1)
        I=find(data.trialStruct.rProfile(ii,19)==model_id);
        model_use_val(ii)=model_val(I);
        model_use_dom(ii)=model_dom(I);
        model_use_cat(ii)=model_cat(I);
        model_use_siz(ii)=model_siz(I);
        model_use_lfq(ii)=model_lfq(I);
        model_use_ctc(ii)=model_ctc(I);
        model_use_trl(ii)=ii;
        model_use_nlt(ii)=model_nlt(I);
    end
    
    % This has resulted in weird dimenisons I think, add some safety
    model_use_trl=zscore(model_use_trl)';
    if size(model_use_trl,2)==1,
        model_use_trl=model_use_trl';
    end
    
    switch model_name
        case 'COCA_v2'
            X=[model_use_val(~reject_val)' ...
                mem_word_val(~reject_val)' ...
                mem_pict_val(~reject_val)' ...
                mem_pict_val(~reject_val)'.*mem_word_val(~reject_val)' ...
                model_use_trl(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_lfq(~reject_val)'.*model_use_val(~reject_val)' ...
                model_use_siz(~reject_val)' ...
                model_use_cat(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=20;
       case 'COCA_v3'
            X=[model_use_val(~reject_val)' ...
                mem_word_val(~reject_val)' ...
                mem_pict_val(~reject_val)' ...
                mem_pict_val(~reject_val)'.*mem_word_val(~reject_val)' ...
                model_use_trl(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_lfq(~reject_val)'.*model_use_val(~reject_val)' ...
                model_use_siz(~reject_val)' ...
                model_use_cat(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=20;
        case 'COCA_v4'
            X=[model_use_val(~reject_val)' ...
                mem_word_val(~reject_val)' ...
                mem_pict_val(~reject_val)' ...
                mem_pict_val(~reject_val)'.*mem_word_val(~reject_val)' ...
                model_use_trl(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_lfq(~reject_val)'.*model_use_val(~reject_val)' ...
                model_use_siz(~reject_val)' ...
                model_use_cat(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=20;
        case 'COCA_v5'
            X=[mem_word_val(~reject_val)' ...
                mem_pict_val(~reject_val)' ...
                mem_pict_val(~reject_val)'.*mem_word_val(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=4;
        case 'COCA_v6'
             X=[mem_word_val(~reject_val)' ...
                model_use_val(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_lfq(~reject_val)'.*model_use_val(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=5;     
        case 'COCA_v7'
             X=[mem_word_val(~reject_val)' ...
                model_use_val(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_lfq(~reject_val)'.*model_use_val(~reject_val)'];
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
            glm_mod.Nv=5; 
        case 'COCA_v8'
%              X=[mem_word_val(~reject_val)' ...
%                 model_use_val(~reject_val)' ...
%                 model_use_lfq(~reject_val)'];
%             glm_mod.Nv=4; 
            
            
             X=[mem_word_val(~reject_val)' ...
                model_use_val(~reject_val)' ...
                model_use_lfq(~reject_val)' ...
                model_use_val(~reject_val)'.*mem_word_val(~reject_val)' ...
                model_use_lfq(~reject_val)'.*mem_word_val(~reject_val)' ...
                model_use_val(~reject_val)'.*model_use_lfq(~reject_val)'];
            glm_mod.Nv=7;
            % Not really using this anymore, setting up myself, under the
            % asssumption that the last input is the categorical var
        case 'COCA_v9'
            X=[mem_word_val(~reject_val)' ...
                model_use_val(~reject_val)' ...
                model_use_lfq(~reject_val)'];
            glm_mod.Nv=4; 
    end

    %=======================================%
    % 3) Run some models
    %=======================================%
    % Smooth the data
    if smooth_data==1
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 10; % Effectively 100ms moving average
        data = ft_preprocessing(cfg,data); 
        % Nothing over 10Hz now anyways, let's downsample to speed up the
        % process...
        if glm_mod.ds==1
            cfg_ds.resamplefs=100;
            cfg_ds.detrend='no';
            cfg_ds.demean='no';
            data=ft_resampledata(cfg_ds,data);
        end
    end
    % Run the base model
    [T,est,timecourse,R]=eeg_glm(X,data.trial(logical(~reject_val)),1:64,251:625,glm_mod);
    save(glm_save,'T','est','timecourse','R','X');
    clear T est timecourse R;
    
    if glm_mod.freq==1
        for ii=1:length(Fband)
            [T,est,timecourse,R]=eeg_glm(X,Fband{ii}(logical(~reject_val)),[],[],glm_mod);
            save(fullfile(RUN.adv.save,[model_name '_' Band_names{ii} '_' RUN.dir.subjects{iSubj} '.mat'])...
                ,'T','est','timecourse','R');
        end
    end
end
end
%=========================================================================%
%% Evaluating Results
%=========================================================================%
model_coef=glm_mod.model_coef;
% data must be loaded here...
% Category 8 appears to be the reference, as it is listed as the first
% trial in the X array, but, this needs closer examined, as if the first
% trial is rejected then this value changes
if glm_mod.freq==1
    LV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
    A=load(LV);
    ga_template.label=A.freqlock.label;
    ga_template.freq=cfg_freq.foi;
    ga_template.time=cfg_freq.toi;
    ga_template.dimord='subj_chan_freq_time';
    clear A;
    % ga_freq.(label).tempalte
    % .powspctrm(ID,ch,freq,time)
end

for iSubj=1:length(RUN.dir.subjects)
    TIME=fullfile(RUN.adv.save,[model_name '_' RUN.dir.subjects{iSubj} '.mat']);
    load(TIME);
    glm_T{iSubj}=T;
    glm_E{iSubj}=est;
    Y{iSubj}=timecourse;
    Rsqr{iSubj}=R;
    design{iSubj}=X;
    clear T est timecourse R X;

    % Need a freq template to steal from
    if glm_mod.freq==1
        for ii=1:length(Band_names)
            FREQ=fullfile(RUN.adv.save,[model_name '_' Band_names{ii} '_' RUN.dir.subjects{iSubj} '.mat']);
            load(FREQ);
            for kk=1:length(model_coef)
                ga_freq.([model_coef{kk} '_v']).label=ga_template.label;
                ga_freq.([model_coef{kk} '_v']).freq=ga_template.freq;
                ga_freq.([model_coef{kk} '_v']).time=ga_template.time;
                ga_freq.([model_coef{kk} '_v']).dimord=ga_template.dimord;
                ga_freq.([model_coef{kk} '_v']).powspctrm(iSubj,:,ii,:)=est{kk};

                ga_freq.([model_coef{kk} '_T']).label=ga_template.label;
                ga_freq.([model_coef{kk} '_T']).freq=ga_template.freq;
                ga_freq.([model_coef{kk} '_T']).time=ga_template.time;
                ga_freq.([model_coef{kk} '_T']).dimord=ga_template.dimord;
                ga_freq.([model_coef{kk} '_T']).powspctrm(iSubj,:,ii,:)=T{kk};

                ga_freq.([model_coef{kk} '_t']).label=ga_template.label;
                ga_freq.([model_coef{kk} '_t']).freq=ga_template.freq;
                ga_freq.([model_coef{kk} '_t']).time=ga_template.time;
                ga_freq.([model_coef{kk} '_t']).dimord=ga_template.dimord;
                ga_freq.([model_coef{kk} '_t']).powspctrm(iSubj,:,ii,:)=timecourse;

                ga_freq.([model_coef{kk} '_R2']).label=ga_template.label;
                ga_freq.([model_coef{kk} '_R2']).freq=ga_template.freq;
                ga_freq.([model_coef{kk} '_R2']).time=ga_template.time;
                ga_freq.([model_coef{kk} '_R2']).dimord=ga_template.dimord;
                ga_freq.([model_coef{kk} '_R2']).powspctrm(iSubj,:,ii,:)=R;
            end
            clear T est timecourse R;
        end
    end
end

LV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
load(LV);
[gaM,gaT,gaS]=eeg_GLM_group(glm_E,glm_T,Y,Rsqr,design,model_coef);

% Setup the templates
sdisp('Plotting',1);
load(RUN.template.plt);
template{2}=layoutmw64; clear layoutmw64;
template{2}.scale=[1 1];
load(RUN.template.plt2);
template{1}=layoutmw64_martyPlot; clear layoutmw64_martyPlot;
template{1}.scale=[2.3 1.5];

RUN.dir.QAL=fullfile(RUN.dir.QAL,'zzRESULTS',RUN.dir.sess,'GLM');
if ~exist(RUN.dir.QAL,'dir'), mkdir(RUN.dir.QAL); end

% Ghetto fix freq information
fn=fieldnames(ga_freq);
for jj=1:length(fn)
    ga_nfreq.(fn{jj})=ga_freq.(fn{jj});
    s=size(ga_freq.(fn{jj}).powspctrm);
    altr=1:2:s(3); 
    ga_nfreq.(fn{jj}).powspctrm=ga_freq.(fn{jj}).powspctrm(:,:,altr,:); 
    ga_nfreq.(fn{jj}).freq=ga_freq.(fn{jj}).freq(1:2:end);   
end
ft_plot_freq(ga_nfreq,'WordM_v','WordM_v',template,0,model_name,'test')
ft_plot_freq(ga_nfreq,'Freq_v','Freq_v',template,0,model_name,'test')
ft_plot_freq(ga_nfreq,'LFreq_v','LFreq_v',template,0,model_name,'test')

ft_plot(gaM,'High_WordM','Low_WordM',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);
ft_plot(gaM,'High_Freq','Low_Freq',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);
% 
A=gaM.High_Freq;
B=gaM.Low_Freq;

A=gaM.High_WordM;
B=gaM.Low_WordM;

figure(1);
subplot(1,4,1);
O1=erp_ss_plot(A,B,[23 27 29],[0.250 0.400],'LF','Hit vs Miss [250ms:400ms]');
subplot(1,4,2);
O2=erp_ss_plot(A,B,[15 21],[0.250 0.400],'LPF','Hit vs Miss [250ms:400ms]');
subplot(1,4,3);
O3=erp_ss_plot(A,B,[24 28 30],[0.250 0.400],'RF','Hit vs Miss [250ms:400ms]');
subplot(1,4,4);
O4=erp_ss_plot(A,B,[16 22],[0.250 0.400],'RPF','Hit vs Miss [250ms:400ms]');

figure(2);
subplot(1,2,1);
O3=erp_ss_plot(A,B,[1 5 6],[0.250 0.400],'PF','Hit vs Miss [250ms:400ms]');
subplot(1,2,2);
O4=erp_ss_plot(A,B,[4],[0.600 0.800],'Cz','Hit vs Miss [600ms:800ms]');

figure(3);
subplot(1,2,1);
O5=erp_ss_plot(A,B,[15 21],[0.250 0.400],'LPF','Hit vs Miss [250ms:400ms]');
subplot(1,2,2);


figure(2);
subplot(1,2,1);
O3=erp_ss_plot(A,B,[1 5 6],[0.800 1],'PF','Hit vs Miss [800ms:1000ms]');
subplot(1,2,2);
O4=erp_ss_plot(A,B,[4],[0.600 0.800],'Cz','Hit vs Miss [600ms:800ms]');
ft_plot(gaM,'High_Freq','Low_Freq',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);


ft_plot_freq(ga_nfreq,'WordM_v','WordM_v',template,0,model_name,'test')
ft_plot_freq(ga_nfreq,'Freq_v','Freq_v',template,0,model_name,'test')


% Blue  = 1
% Red   = 2
% Green = 3
% Black = 4
% if glm_mod.mem==0
%     % Need to fix weird Mem regressor
%     ga_freq.WordM_v.powspctrm=-ga_freq.WordM_v.powspctrm;
%     gaM.WordM.avg=-gaM.WordM.avg;
%     gaM.WordM.individual=-gaM.WordM.individual;
% end
ft_plot(gaM,'High_WordM','Low_WordM',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);
ft_plot(gaM,'High_Freq','Low_Freq',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);

ft_plot(gaM,'High_Freq','High_Freq',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);
ft_plot(gaM,'Low_Freq','Low_Freq',0,8,[],template,u_save_str,0,[-15 15],model_name,[],[],0);


ft_plot_freq(ga_freq,'WordM_v','WordM_v',template,0,model_name,'test')
ft_plot_freq(ga_freq,'Freq_v','Freq_v',template,0,model_name,'test')
ft_plot_freq(ga_freq,'LFreq_v','LFreq_v',template,0,model_name,'test')
ft_plot_freq(ga_freq,'Freqaxn_v','Freqaxn_v',template,0,model_name,'test')


ft_plot(gaM,'WordM','WordM',0,8,[],template,u_save_str,0,[-15 15],model_name,'Freq','LFreq',0);
ft_plot(gaM,'Y','WordM',0,8,[],template,u_save_str,0,[-15 15],model_name,'Freq','LFreq',0);

ft_plot(gaT,'WordM','WordM',0,8,[],template,u_save_str,0,[-5 5],model_name,[],[],0);


ft_plot(gaM,'Y','Int',0,8,[],template,u_save_str,0,[-10 10],model_name,[],[],0);
ft_plot(gaM,'Y','WordM',0,8,[],template,u_save_str,0,[-10 10],model_name,'PictM','Maxn',0);
ft_plot(gaM,'Y','Freq',0,8,[],template,u_save_str,0,[-10 10],model_name,'LFreq','Freqaxn',0);


ft_plot(gaT,'WordM','PictM',0,8,[],template,u_save_str,0,[-5 5],model_name,'Maxn',[],0);
ft_plot(gaT,'Freq','LFreq',0,8,[],template,u_save_str,0,[-5 5],model_name,'Freqaxn',[],0);
ft_plot(gaM,'WordM_Rsqr','PictM_Rsqr',0,8,[],template,u_save_str,0,[0 0.2],model_name,'Maxn_Rsqr',[],0);


ft_plot(gaT,'Maxn','Freq',0,8,[],template,u_save_str,1,[-5 5],model_name,'LFreq','Freqaxn',0);
ft_plot(gaM,'Y','WordM',0,8,[],template,u_save_str,0,[-10 10],model_name,'PictM','Maxn',0);
ft_plot(gaM,'R','R',0,8,[],template,u_save_str,0,[0 0.2],model_name,[],[],0);
ft_plot(gaT,'Maxn','Maxn',0,8,[],template,u_save_str,3,[-5 5],model_name,[],[],0);



ft_plot(gaT,'ctg_4','ctg_1',0,8,[],template,u_save_str,0,[-5 5],model_coef{ii},'ctg_2','ctg_3',0);

ft_plot(B3_gaT,'Trial','WordM',0,8,[],template,u_save_str,0,[-5 5],model_coef{ii},'PictM','Maxn',0);


iSubj=0;
for ii=2:length(model_coef)
    ft_plot(gaM,model_coef{ii},'Y',0,8,[],template,u_save_str,iSubj,[-12 12],[model_name '_' model_coef{ii}],[],[],0);
end

for jj=1:length(Band_names)
            ft_plot(eval([Band_names{jj} '_gaM']),'Y',model_coef{ii},0,8,[],...
                template,u_save_str,iSubj,[-1 1],[Band_names{jj} '_' model_coef{ii}],[],[],1);
end


ft_plot(gaM,'Int','PictM',0,8,[],template,u_save_str,0,[-15 15],model_name,[])
ft_plot(gaM,'Freq','WordM',0,8,[],template,u_save_str,0,[-12 12],model_name,'PictM','Maxn')
ft_plot(gaM,'Freq','WordM',0,8,[],template,u_save_str,0,[-12 12],model_name,'PictM','Int')

% model_coef{end+1}='TimeCourse';
% % Ghetto plotting
% x=-400:4:2000; ch=2; cmat={'k:','b','g','r','c','m','k--'};
% for ii=1:5
%     plot(x,T{ii}(ch,:),cmat{ii},'linewidth',4); hold on;
% end
% plot(x,GYv(ch,:),cmat{end},'linewidth',4)
% legend(model_coef); grid;
% 
% plot(x,glm_T{8}{1}','linewidth',4); grid;
% legend(model_coef);


end
%=========================================================================%
%% Coherency
%=========================================================================%
if RUN.adv.coh==1

save_dir=fullfile(RUN.dir.QAL,'zzResults','Coherence',RUN.dir.sess);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

load(gcohs_save); % Loads in ga_coh
fn=fieldnames(ga_coh);
for fwant=1:3
    for twant=1:18
        for f=1:length(fn)
            save_str=['theta_f' num2str(fwant) '_t' num2str(twant)];
            d=ga_coh.(fn{f});
            if ~isempty(d)
                if ~exist(fullfile(save_dir,fn{f}),'dir'),
                    mkdir(fullfile(save_dir,fn{f}));
                end
                S=d.V(:,fwant,twant);
                for s=1:length(S),
                    if s==1, 
                        A=S{s};
                    else
                        A=A+S{s};
                    end
                end
                Av{f}=A./length(S);
                RSA_pdm_GT(Av{f},1,'',...
                    fullfile(save_dir,fn{f},[save_str '.png']));
            end
        end
        % Make some simple contrast
        RSA_pdm_GT(Av{12}-Av{13},0.25,'',...
            fullfile(save_dir,'R1hit_vs_miss',[save_str '.png']));
        RSA_pdm_GT(Av{3}-Av{4},0.25,'',...
            fullfile(save_dir,'R2hit_vs_miss',[save_str '.png']));
        RSA_pdm_GT(Av{9}-Av{10},0.25,'',...
            fullfile(save_dir,'OO_vs_NN',[save_str '.png']));
        RSA_pdm_GT(Av{7}-Av{8},0.25,'',...
            fullfile(save_dir,'ON_vs_NO',[save_str '.png']));
        clear Av A;
    end
end


end % If statement to enter
%=========================================================================%
%% RSA
%=========================================================================%
% for ii=1:length(RUN.dir.subjects)
%     [~,save_prefix,~]=fileparts(['subj_' RUN.dir.inc_subjects{iSubj} '_' RUN.save_str{5}]); 
%     data_file=fullfile(RUN.dir.pro,RUN.dir.inc_subjects{iSubj},pdir,[save_prefix '_timelock.mat']);
%     data_file_freq=fullfile(RUN.dir.pro,RUN.dir.inc_subjects{iSubj},pdir,[save_prefix '_freqlock_' F_method '.mat']);  
% 
%     load(data_file);
%     load(data_file_freq);
%     
%     % Only 1.1 seconds
%     X=timelock.pre_p_1_stm{1}.trial(:,:,1:400);
%     S=size(X);
%     Z=reshape(X,[size(X,1),size(X,2)*size(X,3)]);
%     C=corr(Z');
%     
%     load(fullfile(RUN.dir.pro,RUN.dir.inc_subjects{iSubj},...
%         pdir,[save_prefix '.mat']),'data');
%     eeg_rsa(data,iSubj);
% end

% if RSA==1
% save_dir='D:\Data\SEFER\EEG\rsa\subject_templates\v1_20140807\standard\';
% subj={'2764' '2772' '2781' '2789'};
% save_prefix='_eegepoch_enc_interp_mastoid_ica_complex_f5_hc_stm_SVM_hr.mat';
% X=(1:625)/250;
% 
% for ii=1:length(subj)
%     load(fullfile(save_dir,subj{ii},['subj_' subj{ii} save_prefix]));
%     for jj=1:length(acc),
%         SET{jj}(ii,:)=mean(acc{jj});
%     end
% end
% 
% S=1;
% cmat=jet(4); figure;
% plot(X,smooth(nanmean(SET{1}),S)','linewidth',3,'color',cmat(1,:)); hold on;
% plot(X,smooth(nanmean(SET{2}),S)','linewidth',3,'color',cmat(2,:)); hold on;
% plot(X,smooth(nanmean(SET{3}),S)','linewidth',3,'color',cmat(3,:)); hold on;
% legend('R1','R2','R1R2');
% 
% sf=250;
% L=625;
% NFFT=2^nextpow2(L);
% P=3;
% for ii=1:4
%     figure(1);
%     plot(X,smooth(SET{P}(ii,:),S)','linewidth',3,'color',cmat(ii,:)); hold on;
%     Y=fft(SET{P}(ii,:),NFFT)/L;
%     f=sf/2*linspace(0,1,NFFT/2+1);
%     figure(2);
%     plot(f,2*abs(Y(1:NFFT/2+1)),'linewidth',3,'color',cmat(ii,:)); hold on;
% end
% legend(subj); 
% switch P
%     case 1, title('R1');
%     case 2, title('R2');
%     case 3, title('R1R2');
% end
% 
% end