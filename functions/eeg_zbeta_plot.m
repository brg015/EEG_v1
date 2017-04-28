function gab=eeg_zbeta_plot(BE,BET,mdl_v,dat,data,it,intc,S)
% Model data       => X
% mdl_v            => template model
% dat              => struct dat
% it               => times included

%=========================================================================%
% Reformat data
%=========================================================================%
% ga.(fn).beta=[chan X time]
% ga.(fn).time=time vector (double)
% ga.(fn).label=
% ga.(fn).dimord='subj_chan_time'
Nch=64; Ntp=sum(it);
mdl_var=mdl_v.CoefficientNames;
if intc==1, mdl_var{1}='intercept'; end
for ii=1:length(mdl_var)
    I=findstr(mdl_var{ii},':');
    if ~isempty(I), mdl_var{ii}(I)='X'; end
end

for ii=1:length(dat)
    m1(ii,:)=dat{ii}.mean;
    m2(ii,:)=dat{ii}.std;
end
m1v=mean(m1); m2v=mean(m2);

gab.('pval').beta=reshape(S,[Nch Ntp]);
gab.('pval').time=data.time{1}(it);
gab.('pval').label=data.label;
    
for ii=1:length(mdl_v.CoefficientNames)
    m=BE(:,ii);
    if intc==1
        if ii>1
            % value*B + intercept
            mv1=BE(:,1).*m2v'+m.*m2v'+m1v';
            mv2=BE(:,1).*m2v'-m.*m2v'+m1v';
        else
            mv1=BE(:,1).*m2v'+m1v';
            mv2=m1v;
        end
    else
        mv1=m.*m2v'+m1v';
        mv2=-m.*m2v'+m1v';  
    end
    gab.([mdl_var{ii} '_pos']).beta=reshape(mv1,[Nch Ntp]);
    gab.([mdl_var{ii} '_pos']).time=data.time{1}(it);
    gab.([mdl_var{ii} '_pos']).label=data.label;
    
    gab.([mdl_var{ii} '_neg']).beta=reshape(mv2,[Nch Ntp]);
    gab.([mdl_var{ii} '_neg']).time=data.time{1}(it);
    gab.([mdl_var{ii} '_neg']).label=data.label;
    
    gab.([mdl_var{ii} '_tstat']).beta=reshape(BET(:,ii),[Nch Ntp]);
    gab.([mdl_var{ii} '_tstat']).time=data.time{1}(it);
    gab.([mdl_var{ii} '_tstat']).label=data.label;
end

end
