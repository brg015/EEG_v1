function granger_preprocess(tdata)
global g;

% 1) Convert tdata to X [n x m x N] = 
% [(channels) variables, observations(time), trials()]
for ii=1:length(tdata.trial)
    X(1,:,ii)=tdata.trial{ii}(g.ch1,g.te);
    X(2,:,ii)=tdata.trial{ii}(g.ch2,g.te);
end

% 2) Standardize Trials in time series
for ii=1:length(tdata.trial)
    X(1,:,ii)=(X(1,:,ii)-mean(X(1,:,ii)))./std(X(1,:,ii));
    X(2,:,ii)=(X(2,:,ii)-mean(X(2,:,ii)))./std(X(2,:,ii));
end

% 3) Subtract the ensemble mean - this should always be taken from all
% of the trials I believe
if isempty(g.ens)
    ens.E1=mean(squeeze(X(1,:,:)),2);
    ens.Es1=std(squeeze(X(1,:,:))');
    ens.E2=mean(squeeze(X(2,:,:)),2);
    ens.Es2=std(squeeze(X(1,:,:))');
    g.ens=ens;
end
X(1,:,:)=squeeze(X(1,:,:))-repmat(g.ens.E1,1,size(X,3));
X(2,:,:)=squeeze(X(2,:,:))-repmat(g.ens.E2,1,size(X,3));
% 3) Normalize
X(1,:,:)=squeeze(X(1,:,:))./repmat(g.ens.Es1',1,size(X,3));
X(2,:,:)=squeeze(X(2,:,:))./repmat(g.ens.Es2',1,size(X,3));

% Fix g.tv wrt g.te
I=g.tv;
I(g.te==0)=[];
X=X(:,I,:);

g.nvars=size(X,1);
g.ntrials=size(X,3);
g.nobs=size(X,2);
g.X=X;