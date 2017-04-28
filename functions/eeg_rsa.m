function eeg_rsa(data,iSubj)
%=========================================================================%
% Create RSA templates
%=========================================================================%
global RUN;

[~,save_prefix,~]=fileparts(['subj_' RUN.dir.inc_subjects{iSubj} '_' RUN.save_str{5}]);
%=========================%
% Rm Kill Channels
%=========================%
rm_ch={'LVEOG','LM','RM'};   
for ii=1:length(rm_ch)
    k_ch(ii)=find(strcmp(rm_ch{ii},data.label));
end
save_ch=setdiff(1:64,k_ch);

sdisp('RSA',1);

save_dir=fullfile(RUN.rsa.sdir,RUN.dir.ver,RUN.rsa.analysis,RUN.dir.inc_subjects{iSubj});
if ~exist(save_dir,'dir'), mkdir_tree(save_dir), end

% Focus only on encoding at first:
rej=~or(data.trialStruct.reject.thresh_indx,data.trialStruct.reject.trend_indx);
idx{1}=find(and(data.trialStruct.event==1,rej)); % pre
idx{2}=find(and(data.trialStruct.event==2,rej)); % stm
for ii=1:length(idx)
    % Collect Trials
    for jj=1:length(idx{ii})
        Time{ii}(jj,:,:)=data.trial{idx{ii}(jj)}(save_ch,:);
    end
end

for ii=1:length(idx)
    
sdisp(['Set ' num2str(ii)],2)
%=========================================================================%    
rProfile = data.trialStruct.rProfile(idx{ii},:);
R1hit=rProfile(:,4);
R1miss=rProfile(:,5);
R2hit=rProfile(:,9);
R2miss=rProfile(:,10);
R1R2hit=rProfile(:,15);
R1R2miss=rProfile(:,16);
Niter=5;

L=1;
Label{L} = R1hit; 
Label{L}(Label{L}==1)=2;
Label{L} = Label{L} + R1miss;
Col{L}=Label{L}>0;
Label{L}(Label{L}==0)=[];
Set{L}=Randi(Niter,[1,length(Label{L})]);
L=L+1;

Label{L} = R2hit; 
Label{L}(Label{L}==1)=2;
Label{L} = Label{L} + R2miss;
Col{L}=Label{L}>0;
Label{L}(Label{L}==0)=[];
Set{L} = zeros(1,length(Label{L}));
Set{L}=Randi(Niter,[1,length(Label{L})]);
L=L+1;

Label{L} = R1R2hit; 
Label{L}(Label{L}==1)=2;
Label{L} = Label{L} + R1R2miss;
Col{L}=Label{L}>0;
Label{L}(Label{L}==0)=[];
Set{L} = zeros(1,length(Label{L}));
Set{L}=Randi(Niter,[1,length(Label{L})]);
L=L+1;
%=========================================================================%
    tL=size(Time{ii},3);  
    for aa=1:L-1
        fprintf(['|- Class ' num2str(aa) '\n']);
        for jj=1:tL % Loop through timepoints   
            Data=squeeze(Time{ii}(Col{aa},:,jj));
            for kk=1:Niter % Iterations
                try
                    Itrain=Set{aa}~=kk;
                    Itest=Set{aa}==kk;
                    % Train
                    % Training: [observations X features]
                    % Group:    [label of observations]
                    SVMStruct=svmtrain(Data(Itrain,:),Label{aa}(Itrain),'method','LS');
                    % Test
                    % SVMStruct: from above
                    % Sample:    [observations X features]
                    % Group:     [predicted value]
                    Guess=svmclassify(SVMStruct,Data(Itest,:));
                    % Sensitivity and Specificity
                    zz=1; contra=2; 
                    TP(zz)=sum((Guess==zz)+(Label{aa}(Itest)==zz)==2);
                    FP(zz)=sum((Guess==zz)+(Label{aa}(Itest)==contra)==2);
                    TN(zz)=sum((Guess==contra)+(Label{aa}(Itest)==contra)==2);
                    FN(zz)=sum((Guess==contra)++(Label{aa}(Itest)==zz)==2);
                    MCC(zz)=((TP(zz)*TN(zz)-FP(zz)*FN(zz))) / ...
                        sqrt((TP(zz)+FP(zz))*(TP(zz)+FN(zz))*(TN(zz)+FP(zz))*(TN(zz)+FN(zz)));
                    acc{aa}(kk,jj)=MCC(zz);
                    clear SVMStruct Guess TP FP TN FN;
                catch err
                    keyboard
                    acc{aa}(kk,jj)=NaN;
                end
            end
            clear Data;   
        end
    end
    switch ii
        case 1, save(fullfile(save_dir,[save_prefix '_pre_SVM_hr.mat']),'acc');
        case 2, save(fullfile(save_dir,[save_prefix '_stm_SVM_hr.mat']),'acc');
    end
    clear Label Set Col
end






