function [booleans,RT] = createConditions_v2(trialStruct,iSubj)
global RUN; 
switch RUN.dir.sess
    case 'enc', pre=RUN.enc.pre; phase=1;   trialStruct.phase=ones(1,length(trialStruct.rejthresh));
    case 'ret', 
        switch RUN.dir.retset
            case 2, pre=RUN.ret.pre; phase=2;
            case 3, pre=RUN.ret.pre; phase=3;
        end
    otherwise, error('unrecognized experimental period'); 
end
booleans=[]; rejnone=0;
sdisp('Currently not processing ind. categ.',1);

%=========================================================================%
% Get RT data
%=========================================================================%
A=excel_reader(fullfile('D:\Data\SEFER\behav\eeg\pre',...
    ['subject' RUN.dir.subjects{iSubj} '.csv']));

ID_eeg=trialStruct.rProfile(:,19);  % Phase Specific EEG trials
ID=cell2num(A{1}.col);              % All trials ~430
rej=~trialStruct.rejthresh';        % Rejection on threshold only*
if strcmp(RUN.dir.sess,'ret')
    rej=(rej & trialStruct.phase==2);
end

for ii=1:length(ID_eeg), I(ii)=find(ID_eeg(ii)==ID); end

CatchResp=cell2num(A{15}.col(I),'Nan');
CatchHit=cell2num(A{16}.col(I),'Nan');

R1hc=cell2num(A{18}.col(I),'Nan');
R1hit=cell2num(A{19}.col(I),'Nan');
R1miss=cell2num(A{20}.col(I),'Nan');
R1fa=cell2num(A{21}.col(I),'Nan');
R1cr=cell2num(A{22}.col(I),'Nan');

R2hc=cell2num(A{24}.col(I),'Nan');
R2hit=cell2num(A{25}.col(I),'Nan');
R2miss=cell2num(A{26}.col(I),'Nan');
R2fa=cell2num(A{27}.col(I),'Nan');
R2cr=cell2num(A{28}.col(I),'Nan');

RT0=cell2num(A{17}.col(I),'Nan');
RT1=cell2num(A{23}.col(I),'Nan'); 
RT2=cell2num(A{29}.col(I),'Nan'); 
RT={RT0,RT1,RT2};

for jj=1:2
    I=[]; str1=[]; str2=[];
    switch jj
        case 1, str1='Conc'; I=R1hc;
                H=R1hit; M=R1miss; F=R1fa; C=R1cr;
        case 2, str1='Perc'; I=R2hc;
                H=R2hit; M=R2miss; F=R2fa; C=R2cr;
    end
    for ii=1:3
        switch ii
            case 1, str2=['HC' '_' str1]; Ivect=(rej & I==1);
            case 2, str2=['LC' '_' str1]; Ivect=(rej & I==0);
            case 3, str2=['' '' str1];   Ivect=rej;
        end
        
        booleans.([str2 '_Hit'])=(H & Ivect);
        booleans.([str2 '_Miss'])=(M & Ivect);
        % Compute FA and CR for retrieval
        if strcmp(RUN.dir.sess,'enc'), continue; end
        
        booleans.([str2 '_FA'])=(F & Ivect);
        booleans.([str2 '_CR'])=(C & Ivect);
    end % Confidence
end % Task
booleans.CatchResp=(rej & CatchResp);
booleans.CatchHit=(rej & CatchHit);
booleans.All=rej;

return;
%=========================================================================%
% Define Booleans
%=========================================================================%
% 
% for aa=1:2
%     varlist={['R' num2str(aa) 'hit'],['R' num2str(aa) 'miss'],...
%         ['R' num2str(aa) 'fa'],['R' num2str(aa) 'cr'],...
%         ['RT' num2str(aa)]};
%     % Some lazy loop cooding
%     for ii=1:length(ID_eeg)
%         I=ID_eeg(ii)==ID;
%         for jj=1:length(varlist)
%             eval([varlist{jj} 'v(ii)=' varlist{jj} '(I);']);
%         end
%         Hit_v(ii)=Hit(I);
%         Miss_v(ii)=Miss(I);
%         clear I;
%     end
%     for ii=1:length(varlist)
%         X=eval([varlist{ii} 'v']);
%         % Find fast/slow
%         I=X>=1; M=median(eval([varlist{end} 'v(I);']));
%         eval([varlist{ii} 'v_fast=and(' varlist{end} 'v>=M,' varlist{ii} 'v==1);']);
%         eval([varlist{ii} 'v_slow=and(' varlist{end} 'v<M,' varlist{ii} 'v==1);']);
%         clear X I M;
%     end
% end
% clear A RT1 RT2;
% %=========================================================================%
% % Modify desccription mat for HC and LC trails
% %=========================================================================%
% if ~isempty(cond)
%     if cond==1        % LC condition
%         T1_vect=~trialStruct.rProfile(:,3); % R1_hc
%         T2_vect=~trialStruct.rProfile(:,8); % R2_hc
%     elseif cond==2    % HC condition
%         T1_vect=trialStruct.rProfile(:,3); % R1_hc
%         T2_vect=trialStruct.rProfile(:,8); % R2_hc
%     end
%     % Now we have a modification matrix
%     for ii=[4,5,6,7]
%        trialStruct.rProfile(:,ii)=and(trialStruct.rProfile(:,ii),T1_vect); 
%     end
%     for ii=[9,10,11,12]
%        trialStruct.rProfile(:,ii)=and(trialStruct.rProfile(:,ii),T2_vect); 
%     end
% end
% % 1) Let's see how many rejections we get
% %=========================================================================%
% % Determine Trial types
% %=========================================================================%
% % 1) Determine what kind of conditions we'd like to examine. note
% % that catch trials have already been removed in the creation of
% % trialStruct.
% switch RUN.LS
%     case '_strict'
%         rej=~or(trialStruct.reject.thresh_indx,trialStruct.reject.trend_indx); 
%     otherwise
%         rej=~trialStruct.rejthresh;
% end
% 
% RUN.QA_ch{iSubj}.perc=sum(rej)/length(rej);
% 
% if rejnone==1, rej=ones(1,length(rej)); end
% 
% if size(rej,2)>1, rej=rej'; end
% tp={'pre' 'stm' 'rsp'};
% for p=phase
%     display(['Bool Phase: ' num2str(p)]);
%     phase_array=(trialStruct.phase==p)';
%     for ii=2 % Only care about stm period
%         display(['  Period: ' num2str(ii)]);
%         type_match=(trialStruct.event==ii)';
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_stm_all=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_stm_all=and(type_match,and(rej,phase_array));']);
%         stm=type_match-trialStruct.rProfile(:,1); % This is catch responses
%         % They serve as a nice variable to force contrast not to include
%         % catch items e.g. high-freq where behavve doesn't auto remove
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_stm=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_stm=and(stm,and(rej,phase_array));']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Resp=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Resp=and(and(and(type_match,trialStruct.rProfile(:,1)),rej),phase_array);']);
%         
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1hit=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1hit=and(and(and(type_match,trialStruct.rProfile(:,4)),rej),phase_array);']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1miss=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1miss=and(and(and(type_match,trialStruct.rProfile(:,5)),rej),phase_array);']);
% 
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1cr=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1cr=and(and(and(type_match,trialStruct.rProfile(:,6)),rej),phase_array);']);
% % 
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1fa=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1fa=and(and(and(type_match,trialStruct.rProfile(:,7)),rej),phase_array);']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2hit=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2hit=and(and(and(type_match,trialStruct.rProfile(:,9)),rej),phase_array);']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2miss=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2miss=and(and(and(type_match,trialStruct.rProfile(:,10)),rej),phase_array);']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2cr=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2cr=and(and(and(type_match,trialStruct.rProfile(:,11)),rej),phase_array);']);
% 
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2fa=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R2fa=and(and(and(type_match,trialStruct.rProfile(:,12)),rej),phase_array);']);
%         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2on=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2on=and(and(and(type_match,trialStruct.rProfile(:,13)),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2no=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2no=and(and(and(type_match,trialStruct.rProfile(:,14)),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2oo=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2oo=and(and(and(type_match,trialStruct.rProfile(:,15)),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2nn=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_R1R2nn=and(and(and(type_match,trialStruct.rProfile(:,16)),rej),phase_array);']);
%      
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_HighFreq=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_HighFreq=and(and(and(stm,trialStruct.rProfile(:,17)),rej),phase_array);']);
% %         
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_LowFreq=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_LowFreq=and(and(and(stm,trialStruct.rProfile(:,18)),rej),phase_array);']);
% %         
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_High_L_Freq=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_High_L_Freq=and(and(and(stm,trialStruct.rProfile(:,20)),rej),phase_array);']);
% %         
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Low_L_Freq=logical(zeros(1,length(trialStruct.event)));']);
%         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Low_L_Freq=and(and(and(stm,trialStruct.rProfile(:,21)),rej),phase_array);']);
%         
%         booleans.stm_p_1_LL_WL=and(booleans.stm_p_1_Low_L_Freq,booleans.stm_p_1_LowFreq);
%         booleans.stm_p_1_LL_WH=and(booleans.stm_p_1_Low_L_Freq,booleans.stm_p_1_HighFreq);
%         booleans.stm_p_1_LH_WL=and(booleans.stm_p_1_High_L_Freq,booleans.stm_p_1_LowFreq);
%         booleans.stm_p_1_LH_WH=and(booleans.stm_p_1_High_L_Freq,booleans.stm_p_1_HighFreq);
%         booleans.stm_p_1_Con=or(booleans.stm_p_1_LL_WL,booleans.stm_p_1_LH_WH);
%         booleans.stm_p_1_Inc=or(booleans.stm_p_1_LH_WL,booleans.stm_p_1_LL_WH);
% 
% 
%         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Birds=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Birds=and(and(and(stm,trialStruct.rProfile(:,23)==1),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Build=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Build=and(and(and(stm,trialStruct.rProfile(:,23)==2),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Cloth=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Cloth=and(and(and(stm,trialStruct.rProfile(:,23)==3),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Food=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Food=and(and(and(stm,trialStruct.rProfile(:,23)==4),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Fruit=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Fruit=and(and(and(stm,trialStruct.rProfile(:,23)==5),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Furn=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Furn=and(and(and(stm,trialStruct.rProfile(:,23)==6),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Mam=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Mam=and(and(and(stm,trialStruct.rProfile(:,23)==7),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Mus=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Mus=and(and(and(stm,trialStruct.rProfile(:,23)==8),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Street=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Street=and(and(and(stm,trialStruct.rProfile(:,23)==9),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Tool=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Tool=and(and(and(stm,trialStruct.rProfile(:,23)==10),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veg=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veg=and(and(and(stm,trialStruct.rProfile(:,23)==11),rej),phase_array);']);
% %         
% %          eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veh=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veh=and(and(and(stm,trialStruct.rProfile(:,23)==12),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veh=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Veh=and(and(and(stm,trialStruct.rProfile(:,23)==12),rej),phase_array);']);
%         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Living=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Living=and(and(and(stm,trialStruct.rProfile(:,24)==1),rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_NonLiving=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_NonLiving=and(and(and(stm,trialStruct.rProfile(:,24)==2),rej),phase_array);']);
% %       
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Compound=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Compound=and(and(trialStruct.rProfile(:,22)>1,rej),phase_array);']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Hit=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Hit=and(Hit_v'',and(rej,phase_array));']);
% %         
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Miss=logical(zeros(1,length(trialStruct.event)));']);
% %         eval(['booleans.' tp{ii} '_p_' num2str(p) '_Miss=and(Miss_v'',and(rej,phase_array));']);
%         
% %         for aa=1:2
% %             varlist={['R' num2str(aa) 'hit'],['R' num2str(aa) 'miss'],...
% %                 ['R' num2str(aa) 'fa'],['R' num2str(aa) 'cr'],...
% %                 ['RT' num2str(aa)]};
% %             for bb=1:length(varlist)
% %                 for cc=1:2
% %                     switch cc
% %                         case 1, spd='fast';
% %                         case 2, spd='slow';
% %                     end    
% %                     eval(['booleans.' tp{ii} '_p_' num2str(p) '_' varlist{bb} '_' spd '=logical(zeros(1,length(trialStruct.event)));']);
% %                     eval(['booleans.' tp{ii} '_p_' num2str(p) '_' varlist{bb} '_' spd '=and(transpose(' varlist{bb} 'v_' spd '),and(rej,phase_array));']);
% %                 end
% %             end
% %         end
%     end
% end



