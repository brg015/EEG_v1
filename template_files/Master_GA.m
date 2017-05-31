RUN.dir.plot=true(1,38);
RUN.dir.subjects={'3446' '3447','3448','3449','3450',...
    '3451','3452','3453','3454','3455',...
    '3456','3457','3458','3480','3483',...
    '3484','3486','3487','3488','3490',...
    '3491','3492','3500','3501','3507',...
    '3509' '3510','3511','3514','3517',...
    '3519','3520','3525','3526','3528',...
    '3540','3545','3547'}; 
RUN.dir.plot(32)=0;
RUN.dir.plot(25)=0;
RUN.dir.plot(23)=0;
RUN.dir.plot(10)=0;
RUN.dir.plot(9)=0;
RUN.dir.plot(6)=0;
RUN.dir.plot(3)=0;
RUN.dir.plot([34,37,38])=0;
wrk_dir=fullfile('F:\SEA\Enc_EEG\');
RUN.dir.pro=fullfile(wrk_dir,'pro','v0');

save_dir='C:\Users\brg13\Desktop\SEA\Writing\'; N=sum(RUN.dir.plot);
RUN.template.plt3='C:\Users\brg13\Desktop\SEA\layoutmw64.mat';
load(RUN.template.plt3);

c=1; for ii=1:length(RUN.dir.subjects)
    if (RUN.dir.plot(ii)==0), continue; end
    cfg_erp.file{c}=fullfile(RUN.dir.pro,[RUN.dir.subjects{ii} '_timelock.mat']);
    cfg_bool.file{c}=fullfile(RUN.dir.pro,[RUN.dir.subjects{ii} '_booleans.mat']);
    cfg_Fz.file{c}=fullfile(RUN.dir.pro,[RUN.dir.subjects{ii} '_Fz_ENS_DPSS_v3.mat']);
    cfg_LP.file{c}=fullfile(RUN.dir.pro,[RUN.dir.subjects{ii} '_dB_ENS_DPSS_v3.mat']);
c=c+1; end
cfg_erp.type='ERP'; 
cfg_Fz.type='FRQ';
cfg_LP.type='FRQ';

% for ii=1:length(cfg_bool.file)
%     X=load(cfg_bool.file{ii});
%     fn=fieldnames(X.booleans);
%     for jj=1:length(fn)
%         N(ii,jj)=sum(X.booleans.(fn{jj}));
%     end
% end

N=sum(RUN.dir.plot);
ga=ft_ga_custom(cfg_erp);
L={'Rem','RemHi','RemLo','For','ForHi','ForLo','Catch','Targets'}; fn=fieldnames(ga);
for ii=1:length(L)
    ga.([L{ii} '_CI'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+8}),ga.(fn{ii}),'erp'); 
end
save(fullfile(save_dir,['GA_' num2str(N) '_ERP.mat']),'ga'); clear ga;

ga=ft_ga_custom(cfg_Fz);
L={'Rem','RemHi','RemLo','For','ForHi','ForLo','Catch','Targets'}; fn=fieldnames(ga);
for ii=1:length(L)
    ga.([L{ii} '_CI'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+8}),ga.(fn{ii}),'power'); 
end
save(fullfile(save_dir,['GA_' num2str(N) '_Z.mat']),'ga'); clear ga;

ga=ft_ga_custom(cfg_LP);
L={'Rem','RemHi','RemLo','For','ForHi','ForLo','Catch','Targets'}; fn=fieldnames(ga);
for ii=1:length(L)
    ga.([L{ii} '_CI'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+8}),ga.(fn{ii}),'power'); 
end
save(fullfile(save_dir,['GA_' num2str(N) '_LP.mat']),'ga'); clear ga;


