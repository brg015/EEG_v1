function eeg_pull_peaks
global RUN; load layoutmw64.mat;
switch RUN.dir.sess
    case 'enc', pdir='1'; pre.baseline=RUN.enc.pre.baseline;
    case 'ret', pdir='2'; pre.baseline=RUN.ret.pre.baseline;
    case 'err', pdir='3'; pre.baseline=RUN.ret.pre.baseline;
    otherwise, pdir='';
end

Freq=0;
 % Time Window
RUN.peak.ind=0;
RUN.peak.model='LPFC';
RUN.peak.time=[.250 .400];
RUN.peak.elec=[15 21 31];
RUN.peak.freq=[3 6];

% RUN.peak.time=[-.5 1];
% RUN.peak.elec=[60];
% RUN.peak.freq=[3 6];

% RUN.peak.time=[.09 .11];
% RUN.peak.elec=[43 55]; 
% RUN.peak.freq=[3 6];

% RUN.peak.time=[.09 .11];
% RUN.peak.elec=[42 54]; 
% RUN.peak.time=[.09 .11];
% RUN.peak.elec=[43 55 42 54]; 

% RUN.peak.time=[-1.0 1.0];
% RUN.peak.elec=[36 37];
% RUN.peak.freq=[8 20];

% RUN.peak.time=[0.8 1.2];
% RUN.peak.elec=[4 13 14 49 50];
% RUN.peak.time=[0.8 0.9];
% RUN.peak.elec=[44];
% RUN.peak.time=[1.1 1.2];
% RUN.peak.elec=[16 22];


% RUN.peak.time=[-0.5 2.25];


zkill=3; % 3stds
%=========================================================================%
% Load in the COCA data
%=========================================================================%
mdata=excel_reader('D:\Data\SEFER\EEG\template\COCA.csv');
% Consildate
ID=cell2num(mdata{1}.col);
Lemma=zscore(log(cell2num(mdata{7}.col)))';
Lif=cell2num(mdata{13}.col);  % cast to numeric
Cat=cell2num(mdata{14}.col); % cast to numeric
Siz=cell2num(mdata{15}.col);
Lfq=zscore(cell2num(mdata{17}.col))';
Ctc=cell2num(mdata{18}.col);
Cmp=cell2num(mdata{21}.col);
T=[ID;Lemma;Lif;Cat;Siz;Lfq;Ctc;Cmp];

c=1;
for iSubj = 1:length(RUN.dir.subjects) 
    sdisp(RUN.dir.subjects{iSubj},2)
    % Load in temporal data
    post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]; % Based on trials
    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,filesep);
    load(fullfile(save_pro_dir,post_save),'data');

    % Load raw subject data & RT info
    A=excel_reader(fullfile('D:\Data\SEFER\behav\eeg\pre',['subject' RUN.dir.subjects{iSubj} '.csv']));
%         RT1=(cell2num(A{23}.col));
%         RT2=(cell2num(A{29}.col));
%     RT1=(cell2num(A{23}.col)-nanmean(cell2num(A{23}.col)))/nanstd(cell2num(A{23}.col));
%     RT2=(cell2num(A{29}.col)-nanmean(cell2num(A{29}.col)))/nanstd(cell2num(A{29}.col));
    RT1=cell2num(A{23}.col).^-1;
    RT1=(RT1-nanmean(RT1))/nanstd(RT1);
    RT2=cell2num(A{29}.col).^-1;
    RT2=(RT2-nanmean(RT2))/nanstd(RT2);
%     RT1=(cell2num(A{23}.col)-nanmean(cell2num(A{23}.col)));
%     RT2=(cell2num(A{29}.col)-nanmean(cell2num(A{29}.col)));
    RTI=cell2num(A{1}.col);
    clear A;
    
    if Freq==1
        F_method='v5';
        data_file_freq=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
        load(data_file_freq);     

        inc_trials=~data.trialStruct.reject.thresh_indx';
    
        if iSubj==1, Freqband=freqlock.freq; end
        cfg_f=[];
        cfg_f.baselinetype='ztransformSession';
        cfg_f.baseline=[-1 2.5]; % DO NOT CHANGE!
        cfg_f.mu=[];
        cfg_f.sigma=[];
        cfg_f.inc_trials=inc_trials;
        % Ztransform of time data
        Fz=ft_freqbaseline(cfg_f,freqlock);           
        cfg_f=[];
        cfg_f.baselinetype='db';
        cfg_f.baseline=[-1 -0.2]; % DO NOT CHANGE!
        cfg_f.mu=[];
        cfg_f.sigma=[];
        cfg_f.inc_trials=inc_trials;
        FdB=ft_freqbaseline(cfg_f,freqlock); 
    end

    % Loop through all trials
    it=and(data.time{1}>=RUN.peak.time(1),data.time{1}<=RUN.peak.time(2));
    ic=zeros(1,length(data.label)); ic(RUN.peak.elec)=1; ic=logical(ic);
    if Freq==1
        itf=and(Fz.time>=RUN.peak.time(1),Fz.time<=RUN.peak.time(2));
        iff=and(Fz.freq>=RUN.peak.freq(1),Fz.freq<=RUN.peak.freq(2));
    end
    % 0 -> 1.4s
%     A1=squeeze(mean(mean(FdB.powspctrm(:,ic,iff,11:25))));
%     [~,I]=min(mean(A1,2));
%     A2=find(iff,1); 
%     fmax(iSubj)=Fz.freq(I+A2);
%     fmaxI{iSubj}=A1;
%     iff=and(Fz.freq>=Fz.freq(I+A2-1),Fz.freq<=Fz.freq(I+A2+1));

    clear A1 A2 I;

    for ii=1:length(data.trial)
        I=data.trialStruct.rProfile(ii,19);
        dat{iSubj}.v(ii)=mean(mean(data.trial{ii}(ic,it)));    
%         dat{iSubj}.v(ii)=max(max(data.trial{ii}(ic,it))); 
        dat{iSubj}.i(ii)=I;
        dat{iSubj}.c(ii)=data.trialStruct.rProfile(ii,3); % Conf
        dat{iSubj}.c2(ii)=data.trialStruct.rProfile(ii,8); % Conf
        dat{iSubj}.h(ii)=data.trialStruct.rProfile(ii,4); % Hit
        dat{iSubj}.m(ii)=data.trialStruct.rProfile(ii,5); % Miss
        dat{iSubj}.RT(ii)=RT1(RTI==I);
        dat{iSubj}.h2(ii)=data.trialStruct.rProfile(ii,9); % Hit
        dat{iSubj}.m2(ii)=data.trialStruct.rProfile(ii,10); % Miss
        dat{iSubj}.RT2(ii)=RT2(RTI==I);
        % Frequency measures
        if Freq==1
            dat{iSubj}.f_base(ii,:)=squeeze(mean(mean(Fz.powspctrm(ii,ic,iff,itf),2),3));
            dat{iSubj}.f_evkd(ii,:)=squeeze(mean(mean(FdB.powspctrm(ii,ic,iff,itf),2),3));        
        end
    end            
    dat{iSubj}.z=zscore(dat{iSubj}.v);
    dat{iSubj}.zrej=abs(dat{iSubj}.z)>zkill;
    % Correct the zscore based upon outliers
    dat{iSubj}.z(~dat{iSubj}.zrej)=zscore(dat{iSubj}.v(~dat{iSubj}.zrej));

    %=====================================================================%
    % Potential pre-processing
    %=====================================================================%
    % Save ish to X and Y
    for ii=1:length(data.trial)
        % 1) Check the ID    
        I=find(dat{iSubj}.i(ii)==T(1,:));
        BAD=[logical(T(7,I)), logical(T(8,I)==2), dat{iSubj}.zrej(ii)];
        if sum(BAD)==0 % use zrej for now
            % 2) Now we setup those measures...
            P=0; P2=0; P3=0; P4=0;
            if RUN.dir.hc==1
                if dat{iSubj}.h(ii)==1, P=1; end
                if dat{iSubj}.m(ii)==1, P=-1; end
                if dat{iSubj}.h2(ii)==1, P2=1; end
                if dat{iSubj}.m2(ii)==1, P2=-1; end
            elseif RUN.dir.hc==0
                if (dat{iSubj}.m(ii)==1 && dat{iSubj}.c(ii)==1), P=-1.5; end
                if (dat{iSubj}.m(ii)==1 && dat{iSubj}.c(ii)==0), P=-0.5; end       
                if (dat{iSubj}.h(ii)==1 && dat{iSubj}.c(ii)==0), P=0.5; end
                if (dat{iSubj}.h(ii)==1 && dat{iSubj}.c(ii)==1), P=1.5; end
                
                if (dat{iSubj}.m2(ii)==1 && dat{iSubj}.c2(ii)==1), P2=-1.5; end
                if (dat{iSubj}.m2(ii)==1 && dat{iSubj}.c2(ii)==0), P2=-0.5; end       
                if (dat{iSubj}.h2(ii)==1 && dat{iSubj}.c2(ii)==0), P2=0.5; end
                if (dat{iSubj}.h2(ii)==1 && dat{iSubj}.c2(ii)==1), P2=1.5; end
            end
            
            if (P==1 && P2==-1), P3=1; end
            if (P==-1 && P2==1), P4=1; end
 
            switch RUN.peak.ind
                case 0
                    X(c,:)=[iSubj, T(2:6,I)', dat{iSubj}.RT(ii), dat{iSubj}.RT2(ii), P P2 P3 P4]; 
                    Y(c)=dat{iSubj}.z(ii);
                    Yv(c)=dat{iSubj}.v(ii);
                    if Freq==1
                        Fb(c,:)=dat{iSubj}.f_base(ii,:);
                        Fe(c,:)=dat{iSubj}.f_evkd(ii,:);
                    end
                    c=c+1;
                case 1
                    X(c,:)=[iSubj, T(2:6,I)', dat{iSubj}.RT(ii), dat{iSubj}.RT2(ii), P P2 P3 P4]; 
                    Y(c,:)=dat{iSubj}.z(ii,:);    
                case 2
                    if ~exist('X','var'), X=[]; end
                    if ~exist('Y','var'), Y=[]; end
                    A=dat{iSubj}.z(ii,:);
                    Y=[Y,A]; 
                    XT=[iSubj, T(2:6,I)', dat{iSubj}.RT(ii), dat{iSubj}.RT2(ii), P P2 P3 P4];
                    XT=repmat(XT,length(A),1);
                    XT(:,13)=1:length(A);
                    X=[X;XT]; clear XT A;
            end
        end
    end
end
%=========================================================================%
% Subjective Effects
%=========================================================================%
% load('D:\Data\SEFER\EEG\template\cdata.mat','cdata');
% X2=nan(size(X,1),2);
% for ii=1:length(RUN.dir.subjects)
%     sI=find(X(:,1)==ii); % Subject match
%     
%     csI=strcmp(RUN.dir.subjects{ii},cdata.subjects);
%     for jj=1:12
%         % FOI/GK
%         sv=cdata.zdat(csI,[jj,jj+12]);
%         cI=find(X(:,4)==jj); % Category indx
%         I=intersect(sI,cI);
%         
%         X2(I,1:2)=repmat(sv,length(I),1);
%     end
% end

%=========================================================================%
% Linear Model
%=========================================================================%
% Interaction measures
% 1st set => FOI, 2nd set => GK
% Load in template data...
load(RUN.template.plt);
template{2}=layoutmw64; clear layoutmw64;
template{2}.scale=[1 1];
load(RUN.template.plt2);
template{1}=layoutmw64_martyPlot; clear layoutmw64_martyPlot;
template{1}.scale=[2.3 1.5];
%=========================================================================%
% X => [
%   1: Subject 
%   2: WordF 
%   3: Living 
%   4: Cat 
%   5: Size 1=large
%   6: LFreq 
%   7: RT1
%   8: RT2
%   9: M1
%   10: M2
%   11: Alta
%   12: Alta
X2=X;
X2(X2(:,3)==2,3)=0; % Artificial == -1
% Fb -> z
% Fe -> dB
%=========================================================================%
keyboard;

Rdir='D:\Data\SEFER\EEG\R';
% Make vars for each scenario of interest...
% 1&2) clear variables of non-interest
X2(:,[4,5,11,12])=[]; 
X2(:,end+1)=X2(:,7)+X2(:,8);
% X2-> 1:ID 2:WF 3:LV 4:LF 5:R1 6:R2 7:M1 8:M2
I1=X2(:,5)==0;    % R1 nan
I2=X2(:,6)==0;    % R2 nan
I3=X2(:,7)<=0;        % M1 miss
I4=X2(:,8)<=0;        % M2 miss
I5=(X2(:,7)==0 | X2(:,7)<-1); % HC miss
I6=(X2(:,8)==0 | X2(:,8)<-1); % HC miss
% 3) create output matrices (won't change...)
Xhead={'ID' 'WF' 'LV' 'LF' 'RT1' 'RT2' 'M1' 'M2' 'M' 'ZuV' 'uV'};
% 3a) Conceptual Memory
Xa=[X2, Y', Yv'];
Xa(I1,:)=[];
% 3b) Perceptual Memory
Xb=[X2, Y', Yv'];
Xb(I2,:)=[];
% 3c) Conceptual RT
Xc=[X2, Y', Yv'];
Xc(I3,:)=[];
% 3d) Perceptual RT
Xd=[X2, Y', Yv'];
Xd(I4,:)=[];
% 3e) Conceptual RT
Xe=[X2, Y', Yv'];
Xe(I5,:)=[];
% 3f) Perceptual RT
Xf=[X2, Y', Yv'];
Xf(I6,:)=[];
% 3g) General Mem
Xg=[X2, Y', Yv'];

write_csv(Xa,Xhead,fullfile(Rdir,[RUN.peak.model '_CM.csv']));
write_csv(Xb,Xhead,fullfile(Rdir,[RUN.peak.model '_PM.csv']));
write_csv(Xc,Xhead,fullfile(Rdir,[RUN.peak.model '_CT.csv']));
write_csv(Xd,Xhead,fullfile(Rdir,[RUN.peak.model '_PT.csv']));
write_csv(Xe,Xhead,fullfile(Rdir,[RUN.peak.model '_mCT.csv']));
write_csv(Xf,Xhead,fullfile(Rdir,[RUN.peak.model '_mPT.csv']));
write_csv(Xg,Xhead,fullfile(Rdir,[RUN.peak.model '_M.csv']));

% X2-> 1:ID 2:WF 3:LV 4:LF 5:R1 6:R2 7:M1 8:M2 9:M 10:z

varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'RT' 'Z(uV)'};
Xin=Xc(:,[1 2 3 4 5]); 
% Xin(:,end+1)=Xin(:,
Yin=Xc(:,10);
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVar',1)

%-------------------------------------------------------------------------%
% Predict Mem
%-------------------------------------------------------------------------%
%---------------%
% Conceptual Mem
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Z(uV)' 'WxL' 'Conceptual Memory'};
Xin=Xvarin(:,[1 2 3 4 9 10]); 
Yin=Xvarin(:,7);
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%---------------%
% Perceptual Mem
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Z(uV)' 'WxL' 'Perceptual Memory'};
Xin=Xvarin(:,[1 2 3 4 9 10]); 
Yin=Xvarin(:,8);
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%---------------%
% Conceptual RT
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Z(uV)' 'WxL' 'Conceptual RT'};
Xin=Xvarin(:,[1 2 3 4 9 10]); Xin([I1 | I3],:)=[];
Yin=Xvarin(:,5);              Yin([I1 | I3],:)=[];
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%---------------%
% Perceptual RT
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Z(uV)' 'WxL' 'Perceptual RT'};
Xin=Xvarin(:,[1 2 3 4 9 10]); Xin([I2 | I4],:)=[];
Yin=Xvarin(:,6);              Yin([I2 | I4],:)=[];
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%-------------------------------------------------------------------------%
% Predict uV
%-------------------------------------------------------------------------%
for ii=1:18
    display(ii)
    I=Xvarin(:,1)==ii;
    varin={'WordFreq' 'Living' 'LettFreq' 'CM' 'Z(uV)'};
    Xin=Xvarin(:,[2 3 4 7]); 
    Yin=Y(:);
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin)
    a1(ii,:)=mdl_v.Coefficients.Estimate;
end

%---------------%
% Conceptual Mem
%---------------%

%---------------%
% Perceptual Mem
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Perceptual Memory' 'WxL' 'Z(uV)'};
Xin=Xvarin(:,[1 2 3 4 8 10]); 
Yin=Y;
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%---------------%
% Conceptual RT
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Conceptual RT' 'WxL' 'Z(uV)'};
Xin=Xvarin(:,[1 2 3 4 5 10]);  Xin([I1 | I3],:)=[];
Yin=Y;                        Yin([I1 | I3])=[];
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%---------------%
% Perceptual RT
%---------------%
varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Perceptual RT' 'WxL' 'Z(uV)'};
Xin=Xvarin(:,[1 2 3 4 6 10]);  Xin([I2 | I4],:)=[];
Yin=Y;                        Yin([I2 | I4])=[];
mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1)
%-------------------------------------------------------------------------%
% Predict uV(freq)
%-------------------------------------------------------------------------%
for ii=4:size(Fb,2)
    %---------------%
    % Conceptual Mem
    %---------------%
    varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'WxL' 'Z(uV)'};
    Xin=Xa(:,[1 2 3 4 9]); 
    Yin=Fb(:,ii);
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1);
    vWL(ii)=mdl_v.Coefficients.tStat(end);
    vLet(ii)=mdl_v.Coefficients.tStat(end-1);
    vLife(ii)=mdl_v.Coefficients.tStat(end-2);
    vWord(ii)=mdl_v.Coefficients.tStat(end-3);
    v6(ii)=mdl_v.coefTest;
    %---------------%
    % Conceptual Mem
    %---------------%
    varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Conceptual Memory' 'WxL' 'Z(uV)'};
    Xin=Xvarin(:,[1 2 3 4 7 9]); 
    Yin=Fb(:,ii);
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1);
    v1(ii)=mdl_v.Coefficients.tStat(end-1);
    %---------------%
    % Perceptual Mem
    %---------------%
    varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Perceptual Memory' 'WxL' 'Z(uV)'};
    Xin=Xvarin(:,[1 2 3 4 8 9]); 
    Yin=Fb(:,ii);
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1);
    v2(ii)=mdl_v.Coefficients.tStat(end-1);
    %---------------%
    % Conceptual RT
    %---------------%
    varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Conceptual RT' 'WxL' 'Z(uV)'};
    Xin=Xvarin(:,[1 2 3 4 5 9]);  Xin([I1 | I3],:)=[];
    Yin=Fb(:,ii);                        Yin([I1 | I3],:)=[];
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1);
    v3(ii)=mdl_v.Coefficients.tStat(end-1);
    %---------------%
    % Perceptual RT
    %---------------%
    varin={'ID' 'WordFreq' 'Living' 'LettFreq' 'Perceptual RT' 'WxL' 'Z(uV)'};
    Xin=Xvarin(:,[1 2 3 4 6 9]);  Xin([I2 | I4],:)=[];
    Yin=Fb(:,ii);                       Yin([I2 | I4],:)=[];
    mdl_v=fitglm(Xin,Yin,'linear','VarNames',varin,'CategoricalVars',1);
    v4(ii)=mdl_v.Coefficients.tStat(end-1);
end

plot(Fz.time(itf),v1,'b'); hold on; 
plot(Fz.time(itf),v2,'r'); hold on;
plot(Fz.time(itf),v3,'g'); hold on;
plot(Fz.time(itf),v4,'k'); hold on;
legend({'R1' 'R2' 'R1RT' 'R2RT'});