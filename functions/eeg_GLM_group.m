function [gaM,gaT,gaS]=eeg_GLM_group(glm_E,glm_T,Y,Rsqr,X,model_coef,freq)

global RUN;
switch RUN.dir.sess
    case 'enc', pdir='1'; pre.baseline=RUN.enc.pre.baseline;
    case 'ret', pdir='2'; pre.baseline=RUN.ret.pre.baseline;
    case 'err', pdir='3'; pre.baseline=RUN.ret.pre.baseline;
    otherwise, pdir='';
end

% Load in a ga file for reference
[~,u_save_str,~]=fileparts(RUN.save_str{5});
gerps_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N8_erps.mat']);
Z=load(gerps_save);

save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{1},pdir,filesep);
post_save=['subj_' RUN.dir.subjects{1} '_' RUN.save_str{5}]; % Based on trials*
load(fullfile(save_pro_dir,post_save),'data'); % var=data;
   
if exist('freq','var')
    data.time{1}=freq;
end
    
% Sort by subject
Vname= Rsqr{1};

% Fixes IF resampled... this is finicky...
if length(data.time{1})~=size(glm_E{1}{1},2)
    % Resample
    % t1=data.time{1}(1);
    % t2=data.time{1}(end);
    
    t1=data.time{1}(251);
    t2=data.time{1}(625);
    
    tf=linspace(t1,t2,size(glm_E{1}{1},2));
    data.time{1}=tf; % With that offset
end

%=========================================================================%
% Add in additional high/low regressors for pretty plots
%=========================================================================%
N=length(model_coef);
temp=model_coef;
for ii=2:(N*2-1),
    if ii<=N,
        model_coef{ii}=['High_' model_coef{ii}]; 
        
    else
        model_coef{ii}=['Low_' temp{ii-N}];
    end
end
clear temp;

for ii=1:length(RUN.dir.subjects),
    design=X{ii};
    for jj=1:N
        % These are the estimates...
        % Want these to be based upon uV values
        %=================================================================%
        % Essential interpretation section
        %=================================================================%
        v=glm_E{ii}{jj};               % Raw value for given condition
        if jj==1
            % ii=1 is the intercept term
            Gglm_E{jj}(ii,:,:)=v; 
        else
            intcpt=glm_E{ii}{1};           % Intercept term for subject
            reg=design(:,jj-1);            % Design regressor for jj
            ureg=unique(reg(~isnan(reg))); % Unique regressor values
            v=glm_E{ii}{jj};               % Raw value for given condition
            if ureg<5
                hgh_off=max(ureg);
                low_off=min(ureg);
                vh=intcpt+hgh_off.*v;
                vl=intcpt+low_off.*v;
                clear ureg;
            else
                mreg=mean(reg);
                sreg=std(reg);
                vh=intcpt+(mreg+2*sreg).*v;
                vl=intcpt+(mreg-2*sreg).*v;
                clear mreg sreg;
            end
            Gglm_E{jj+N}(ii,:,:)=vl;  % Low
            Gglm_E{jj}(ii,:,:)=vh;    % High
            clear reg ureg vh vl v intcpt;
        end
        % T-values, don't need to play with these at all
        % Same can be said of R2 values
        Gglm_T{jj}(ii,:,:)=glm_T{ii}{jj};
        Gglm_T{jj+N}(ii,:,:)=glm_T{ii}{jj};
        if iscell(Vname), 
            GR{jj}(ii,:,:)=Rsqr{ii}{jj}; 
            GR{jj+N}(ii,:,:)=Rsqr{ii}{jj}; 
        end
    end
    
    GY(ii,:,:)=Y{ii};
    if ~iscell(Vname), 
        GR(ii,:,:)=Rsqr{ii}; 
    end
    
end
GYv=squeeze(mean(GY));
if ~iscell(Vname), GRv=squeeze(mean(GR)); end

%=========================================================================%
% Grand average
%=========================================================================%
for ii=1:length(Gglm_E)
    M{ii}=squeeze(mean(Gglm_E{ii}));
    S{ii}=squeeze(std(Gglm_E{ii}));
    T{ii}=M{ii}./(S{ii}/sqrt(length(RUN.dir.subjects)));    
    if iscell(Vname), GRv{ii}=squeeze(mean(GR{ii})); end
end
%=========================================================================%
% Finish filling in fields
%=========================================================================%
for ii=1:length(model_coef)

    gaM.(model_coef{ii}).avg=M{ii};
    gaM.(model_coef{ii}).individual=Gglm_E{ii};

    gaM.(model_coef{ii}).label=data.label;
    gaM.(model_coef{ii}).time=data.time{1};
    gaM.(model_coef{ii}).dimord='chan_time';
    
    gaT.(model_coef{ii}).avg=T{ii};
    gaT.(model_coef{ii}).individual=Gglm_T{ii};
    gaT.(model_coef{ii}).label=data.label;
    gaT.(model_coef{ii}).time=data.time{1};
    gaT.(model_coef{ii}).dimord='chan_time';
    
    gaS.(model_coef{ii}).avg=S{ii};
    gaS.(model_coef{ii}).label=data.label;
    gaS.(model_coef{ii}).time=data.time{1};
    gaS.(model_coef{ii}).dimord='chan_time';
    
    if iscell(Vname),
        gaM.([model_coef{ii} '_Rsqr']).avg=GRv{ii};
        gaM.([model_coef{ii} '_Rsqr']).label=data.label;
        gaM.([model_coef{ii} '_Rsqr']).time=data.time{1};
        gaM.([model_coef{ii} '_Rsqr']).dimord='chan_time';
    end
end

gaM.Y.avg=GYv;
gaM.Y.individual=GY;
gaM.Y.label=data.label;
gaM.Y.time=data.time{1};
gaM.Y.dimord='chan_time';

if ~iscell(Vname), 
    gaM.R.avg=GRv;
    gaM.R.individual=GR;
    gaM.R.label=data.label;
    gaM.R.time=data.time{1};
    gaM.R.dimord='chan_time';
end