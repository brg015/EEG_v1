function [glm_T,glm_E,Y,R]=eeg_glm(X,trial_data,ch,tp,glm_mod)
% Inputs
%   X: Design Matrix
%   trial_data{trial}(ch X time)
%

% 1) For each time, channel, run the GLM
if (isempty(ch) && isempty(tp))
    s=size(trial_data{1});
    for ii=1:length(trial_data),
        D(ii,:)=reshape(trial_data{ii},1,[]);
    end
else
    s=[length(ch),length(tp)];
    for ii=1:length(trial_data),
        D(ii,:)=reshape(trial_data{ii}(ch,tp),1,[]);
    end
end



% mdl=fitglm(X,D(:,ii));
% Cnames=mdl.CoefficientNames;
% if size(X,2)>3
% C=X(:,end);    X(:,end)=[];
% Z=dummyvar(C);
% vX=[X,Z(:,1:end-1)];
% X=vX;
% end

% 2) Try those GLMs
for ii=1:size(D,2)
    if rem(ii,1000)==0,
        display([' Event ' num2str(ii) ':Analyzed']);
    end

    y=D(:,ii);
    Y(ii)=mean(y);
    if ~isnan(Y(ii))
        if glm_mod.model_ind==0
            mdl=fitglm(X,y); 
            E(ii,1:glm_mod.Nv)=mdl.Coefficients.Estimate;
            T(ii,1:glm_mod.Nv)=mdl.Coefficients.tStat;
            R(ii)=mdl.Rsquared.Ordinary;
            clear mdl y IB;
        else
            for jj=1:size(X,2)
                mdl=fitglm(X(:,jj),y);
                E(ii,jj)=mdl.Coefficients.Estimate(2);
                T(ii,jj)=mdl.Coefficients.tStat(2);
                R(ii,jj)=mdl.Rsquared.Ordinary;
            end
        end
    else
        E(ii,1:glm_mod.Nv)=NaN;
        T(ii,1:glm_mod.Nv)=NaN;
        if glm_mod.model_ind==0
            R(ii)=NaN;
        else
            R(ii,:)=NaN;
        end
    end
end

% 3) Reshape the estimates
% E [Obs X Est]
% T [Obs X Est]
for ii=1:size(E,2),
    glm_T{ii}=reshape(T(:,ii),s);
    glm_E{ii}=reshape(E(:,ii),s);
end
Y=reshape(Y,s);

if glm_mod.model_ind==0
    R=reshape(R,s);
else
    for ii=1:size(R,2),
        Ri{ii}=reshape(R(:,ii),s);
    end
    clear R; R=Ri;
end
    


    
    