function NET=ft_net_analysis(cfg)

GTchan={'ClustCoef','Degree','ecc','Ci'};
GTglob={'lamda','Eglob','radius','diameter','aClustCoef','aDeg','Q','T'};

K=cfg.K; % Electrodes to remove from the network
for ii=1:length(cfg.files)
    % Early_NET2, Late_NET2
    A=load(cfg.files{ii});
    if ii==1,
        NET.time=A.NET.(cfg.fn{1}).time;
        NET.label=A.NET.(cfg.fn{1}).label;
        NET.label(K)=[];
        NET.cfg=cfg;
    end
    
    for jj=1:length(cfg.fn)
        for tt=1:length(NET.time)
%-------------------------------------------------------------------------%
R=A.NET.(cfg.fn{jj}).(cfg.measure); % Network
% Any stat fix
% RT=A.NET.(cfg.fn{jj}).([cfg.measure])./A.NET.(cfg.fn{jj}).([cfg.measure 'sem']);
% Rp=1-tcdf(RT,sum(A.NET.(cfg.fn{jj}).cfg.trials)-1);
% R(Rp<0.001)=0;

if cfg.cont==1, R=squeeze(R(:,:,1,tt)); end
    
R(logical(eye(64)))=0;              % Zero the diag
R(K,:)=[]; R(:,K)=[];               % Remove K channels
R=abs(R);                           % Make positive

% Modularity calculation
for kk=1:10, [Ciz(kk,:), Qz(kk)]=community_louvain(R); end

[Q,I]=max(Qz); Ci=Ciz(I,:); clear Ciz Qz I;
% Reliabe measures from Kuntzelman 2017
C=clustering_coef_wu(R);
D=sum(R); % Ignore diag
L=distance_wei(1./R); L(logical(eye(size(R))))=0;
[lambda,Eglob,ecc,radius,diameter] = charpath(L);
T=transitivity_wu(R);
aD=mean(D);
aC=mean(C);
%-------------------------------------------------------------------------%
% Organize NET output
%-------------------------------------------------------------------------%
NET.glob_measure.lamda.(cfg.fn{jj}).individual(ii,tt)=lambda;
NET.glob_measure.Eglob.(cfg.fn{jj}).individual(ii,tt)=Eglob;
NET.glob_measure.radius.(cfg.fn{jj}).individual(ii,tt)=radius;
NET.glob_measure.diameter.(cfg.fn{jj}).individual(ii,tt)=diameter;
NET.glob_measure.aClustCoef.(cfg.fn{jj}).individual(ii,tt)=aC;
NET.glob_measure.aDeg.(cfg.fn{jj}).individual(ii,tt)=aD;
NET.glob_measure.Q.(cfg.fn{jj}).individual(ii,tt)=Q;
NET.glob_measure.T.(cfg.fn{jj}).individual(ii,tt)=T;

NET.chan_measure.ClustCoef.(cfg.fn{jj}).individual(ii,:,tt)=C;
NET.chan_measure.Degree.(cfg.fn{jj}).individual(ii,:,tt)=D;
NET.chan_measure.ecc.(cfg.fn{jj}).individual(ii,:,tt)=ecc;
NET.chan_measure.Ci.(cfg.fn{jj}).individual(ii,:,tt)=Ci;

NET.network.(cfg.fn{jj})(ii,tt,:,:)=R;
NET.network.dimord='subj_time_chan_chan';
NET.network.time=NET.time;

if ii==1
    for aa=1:length(GTchan)
        NET.chan_measure.(GTchan{aa}).(cfg.fn{jj}).time=NET.time;
        NET.chan_measure.(GTchan{aa}).(cfg.fn{jj}).label=NET.label;
        NET.chan_measure.(GTchan{aa}).(cfg.fn{jj}).dimord='subj_chan_time';
    end
    for aa=1:length(GTglob)
        NET.glob_measure.(GTglob{aa}).(cfg.fn{jj}).time=NET.time;
        NET.glob_measure.(GTglob{aa}).(cfg.fn{jj}).label=NET.label;
        NET.glob_measure.(GTglob{aa}).(cfg.fn{jj}).dimord='subj_time'; 
    end
end
%-------------------------------------------------------------------------%
clear Q Ci C D L lambda efficiency ecc radius diamter T aD aC;
%-------------------------------------------------------------------------%
        end % Time
    end % Contrast
end % Subject

