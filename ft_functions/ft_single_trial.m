function ft_single_trial(ST)
% ST.subjects           -> list of subjects
% ST.dat{ii}.type       -> ERP or OSC
%           .mir        -> mirror electrodes
%           .time       -> time to pull
%           .elec       -> electrodes to pull
% ST.dat{contrast}.file{subject}
% ST.dat{contrast}.beh{subject}
% ST.ERS                -> ==1 if ERS is on
%-------------------------------------------------------------------------%
% Extract data
%-------------------------------------------------------------------------%
for jj=1:length(ST.subjects)
display(['Running: ' ST.subjects{jj}]);
for ii=1:length(ST.dat)
%-------------------------------------------------------------------------%   
switch ST.dat{ii}.type
    case 'ERP'
        if isempty(rdata)
            Hdata=load(ST.dat{ii}.file{jj}); % ERP data  
            % this takes some time, so save rdata
            for kk=1:length(Hdata.data.trial), rdata(kk,:,:)=Hdata.data.trial{kk}; end
        end
        data=ft_pull_ERP(Hdata.data,ST.dat{ii}.elec,ST.dat{ii}.time,rdata);
        if isfield(ST.dat{ii},'mir')
            data2=ft_pull_ERP(Hdata.data,ST.dat{ii}.mir,ST.dat{ii}.time,rdata);
        end
    case 'OSC'
        if isempty(fdata), fdata=load(ST.dat{2}.file{jj}); end % FRQ data
        data=ft_pull_FRQ(fdata,ST.dat{ii}.elec,ST.dat{ii}.time,ST.dat{ii}.freq);
        if isfield(ST.dat{ii},'mir')
            data2=ft_pull_FRQ(fdata,ST.dat{ii}.mir,ST.dat{ii}.time,ST.dat{ii}.freq);
        end
end
% jj is subject
% ii is contrast

% at this point
% data is a vector of extracted values
% data2 is a vector of extracted values on the opposite side
%-------------------------------------------------------------------------%
% Mirror data if needed
%-------------------------------------------------------------------------%
if isfield(ST.dat{ii},'mir')
    % Contra minus Ipsi extraction - dependent upon data being
    % defined via left electrodes (contra convention)
    %
    % a1 left-electrodes w/ contra targets - right-electrodes w/
    % ipsi targets
    % a2 right-electrodes w/ contra targets - left-electrodes w/
    % ispi targets
    a1=data.*booleans.Right_Targets'-data2.*booleans.Right_Targets';
    a2=data2.*booleans.Left_Targets'-data.*booleans.Left_Targets';
    E{ii}.val=(a1+a2)';
else
    E{ii}.val=data';
end

%-------------------------------------------------------------------------%
% Load behave
%-------------------------------------------------------------------------%
load(ST.dat{ii}.beh{jj});

% Assumes has RT and targets var here
if exist('RT','var'), E{ii}.RT=RT; end
if ERS==1
    % fix faulty var name that crops up
    if ~exist('targets','var'), targets=target; clear target; end
    E{ii}.targets=targets;
end

%----------------------%
% Code response outpt
%----------------------%
E{ii}.trialN=zscore(1:length(booleans.Rem));
E{ii}.Mem=zeros(1,length(booleans.Rem));
E{ii}.Mem(booleans.RemHi)=1.5;
E{ii}.Mem(booleans.RemLo)=0.5;
E{ii}.Mem(booleans.ForLo)=-0.5;
E{ii}.Mem(booleans.ForHi)=-1.5;
if ST.dat{ii}.E==2
    R=booleans.Right_Rem+booleans.Right_For;
    L=booleans.Left_Rem+booleans.Left_For;
    booleans.Right_Targets=logical(R);
    booleans.Left_Targets=logical(L);
end
clear L R;

K=E{ii}.Mem==0;
E{ii}.val(K)=[];
E{ii}.RT(K)=[];
E{ii}.targets(K)=[];
E{ii}.Mem(K)=[];      

E{ii}.val=E{ii}.val;
E{ii}.zval=zscore(E{ii}.val);

E{ii}.R=zeros(1,length(booleans.Right_Targets)); E{ii}.R(booleans.Right_Targets)=1;
E{ii}.L=zeros(1,length(booleans.Left_Targets)); E{ii}.L(booleans.Left_Targets)=1;
E{ii}.R(K)=[];
E{ii}.L(K)=[];
E{ii}.trialN(K)=[];
clear K;
%=========================================================================%
end % data loop
%=========================================================================%
%-------------------------------------------------------------------------%
% If ERS
%-------------------------------------------------------------------------%
if ST.ERS==1
    % Setup a basic sorting to ensure things align
    Eind=E{1}.targets;   % assume start is E
    Rind=E{end}.targets; % and that the end is R
    [C,IE,IR]=intersect(Eind,Rind);
    fn=fieldnames(E{1}); % any index should work here
    for kk=1:length(E)
        if ST.dat{kk}.E==1, I=IE; else I=IR; end
        for aa=1:length(fn)
            E{kk}.(fn{aa})=E{kk}.(fn{aa})(I);
        end
    end
    % ^ this forces them to match...
end
    %=====================================================================%
    % Save data to template
    %=====================================================================%
    E{1}.L(E{1}.L==0)=-1;
    E{1}.R(E{1}.R==0)=-1;
    
    sdata{1}.header='ID';
    for kk=1:length(E{1}.val), sdata{1}.col{kk}=ST.subjects{jj}; end 
    sdata{2}.header='Z_invEncRT';
    sdata{2}.col=zscore(1./E{1}.RT);
    sdata{3}.header='Z_invRetRT';
    sdata{3}.col=zscore(1./E{end}.RT);
    sdata{4}.header='targets';
    sdata{4}.col=E{1}.targets; % any one would work
    sdata{5}.header='Mem';
    sdata{5}.col=E{1}.Mem;
    sdata{6}.header='Right';
    sdata{6}.col=E{1}.R;
    sdata{7}.header='Left';
    sdata{7}.col=E{1}.L;
    sdata{8}.header='Etrial';
    sdata{8}.col=E{1}.trialN;
    sdata{9}.header='Rtrial';
    sdata{9}.col=E{end}.trialN;
    c=1;
    for kk=1:length(E)
        sdata{7+c}.header=ST.dat{kk}.name;
        sdata{7+c}.col=E{kk}.val;
        sdata{7+c+1}.header=['Z_' ST.dat{kk}.name];
        sdata{7+c+1}.col=E{kk}.zval;
        c=c+2;
    end
    %=====================================================================%
    % Clean template by removing outliers and non-responses
    %=====================================================================%
    Fdata{jj}=sdata;
    
     % Clean up looped variables
    clear E sdata targets;
end % subject loop

% Reformat Fdata for saving
for ii=1:length(Fdata)
    if ii==1
        out=Fdata{1};
    else
        for jj=1:length(Fdata{ii})
            out{jj}.col=[out{jj}.col, Fdata{ii}{jj}.col];
        end
    end
end

N=length(out);
out{N+1}.header='binMem';
out{N+1}.col=out{5}.col;
out{N+1}.col(out{N+1}.col>0)=1;
out{N+1}.col(out{N+1}.col<0)=-1;

write_struct(out,'F:\Data2\SEA\R\N25_zscore_mem_022117_dpss.csv');