set_study
%Above line clears all variables

dataPath ='C:\Users\adr25\Desktop\data';
eeglabDir = 'matlab_apps\eeglab13_2_2b\';
addpath('matlab_apps/functions')
addpath('matlab_apps\fieldtrip-20140618\')
addpath('matlab_apps\plot2svg_20120915\')

%what does this do?
rand('state',2000);

channels={
    '5' '6';
    '7' '8';
    '9' '10';
    '11' '12';
    '13' '14';
    '15' '16';
    '17' '18';
    '19' '20';
    '21' '22';
    '23' '24';
    '25' '26';
    '27' '28';
    '29' '30';
    %33 is left mastoid
    '39' '40';
    '41' '42';
    '43' '44';
    '45' '46';
    '47' '48';
    '49' '50';
    '51' '52';
    '53' '54';
    '55' '56';
    '57' '58';
    '59' '60';
    '61' '62';
    '63' '64';
}



%% ERP analysis 
for j=1:length(subjectID);
    pp = subjectID{j};
     
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['participant: ' int2str(j)])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')

    load(fullfile(dataPath,pp, [pp '.mat']));
    
    booleans = createConditions(data.trialStruct);
    fn = fieldnames(booleans);
    
    for iConditions = 1 : length(fn)
        cfg=[];
        cfg.trials = booleans.(fn{iConditions});
        timelock.(fn{iConditions}){j} = ft_timelockanalysis (cfg, data);

    end
end

%%
fn = fieldnames(timelock);
for iConditions = 1 : length(fn)
    
    for iSubj = 1:length(timelock.(fn{iConditions}))
        cfg = [];
        %apply filter
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 30;
        temp.(fn{iConditions}){iSubj} = ft_preprocessing(cfg, timelock.(fn{iConditions}){iSubj});
    end
    
    
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.baseline = [-0.2 0];
             
    %grand average
    ga.(fn{iConditions}) = ft_timelockgrandaverage(cfg, temp.(fn{iConditions}){:});
   

end

%%
%create contrasts

ga.congruent = ga.cueLeft;
ga.congruent.individual = 1/4*(ga.cong_v_L.individual + ga.cong_v_R.individual + ...
    ga.cong_i_L.individual + ga.cong_i_R.individual);

ga.incongruent = ga.cueLeft;
ga.incongruent.individual = 1/4 *(ga.incong_v_L.individual + ga.incong_v_R.individual + ...
    ga.incong_i_L.individual + ga.incong_i_R.individual);
    

contrasts.incongruency = ga.cueLeft;
contrasts.incongruency.individual = ga.incongruent.individual - ga.congruent.individual;


ga.valid = ga.cueLeft
ga.valid.individual = 1/4*(ga.cong_v_L.individual + ga.cong_v_R.individual + ...
    ga.incong_v_L.individual + ga.incong_v_R.individual);

ga.invalid = ga.cueLeft;
ga.invalid.individual = 1/4*(ga.cong_i_L.individual + ga.cong_i_R.individual + ...
    ga.incong_i_L.individual + ga.incong_i_R.individual);


contrasts.validity = ga.cueLeft;
contrasts.validity.individual = ga.invalid.individual - ga.valid.individual;

%left and right
%rexpand this contrast
% contrasts.attend=ga.cueLeft;
% contrasts.attend.individual = ga.cong_v.individual - ga.cong_i.individual;

%example: [grandavg] = ft_timelockgrandaverage(cfg, avg1, avg2, avg3, ...)
ga.cong_v = ga.cueLeft;
ga.cong_v.individual = 1/2*(ga.cong_v_L.individual+  ga.cong_v_R.individual);
ga.incong_v = ga.cueLeft;
ga.incong_v.individual = 1/2*(ga.incong_v_L.individual+  ga.incong_v_R.individual);

ga.cong_i = ga.cueLeft;
ga.cong_i.individual =  1/2*(ga.cong_i_L.individual+  ga.cong_i_R.individual);
ga.incong_i = ga.cueLeft;
ga.incong_i.individual =  1/2*(ga.incong_i_L.individual+  ga.incong_i_R.individual);


%validity effects for congruent trials
contrasts.validity_C = ga.cueLeft;
contrasts.validity_C.individual = ga.cong_i.individual - ...
    ga.cong_v.individual;

contrasts.validity_I = ga.cueLeft;
contrasts.validity_I.individual = ga.incong_i.individual - ...
    ga.incong_v.individual;


%Incongruency Effect for valid trials
contrasts.incongruency_V = ga.cueLeft;
contrasts.incongruency_V.individual = ga.incong_v.individual - ...
    ga.cong_v.individual;

%Incongruency Effect for invalid trials
contrasts.incongruency_I = ga.cueLeft;
contrasts.incongruency_I.individual = ga.incong_i.individual - ga.cong_i.individual;


%plot difference wave for each above! 

%%
%CONTRA - IPSI

%example: ga.conMinIpsiFast = contraMinusIpsi(electrodes,ga.cueLeftFast,ga.rightFast, 'erp');
ga.conMinIpsiCue = contraMinusIpsi(electrodes, ga.cueLeft, ga.cueRight, 'erp');

ga.conMinIpsiCongruentValid = contraMinusIpsi(electrodes, ga.cong_v_L, ga.cong_V_R, 'erp');
ga.conMinIpsiCongruentInvalid = contraMinusIpsi(electrodes, ga.cong_i_L, ga.cong_i_R, 'erp');

ga.conMinIpsiIncongruentValid = contraMinusIpsi(electrodes, ga.incong_v_L, ga.incong_v_R, 'erp');
ga.conMinIpsiIncongruentInvalid = contraMinusIpsi(electrodes, ga.incong_i_L, ga.incong_i_R, 'erp');

contrasts.conMinIpsi = ga.cueLeft;
%double check below lines:
contrasts.conMinIpsi.individual = ga.conMinIpsiCue.individual;

contrasts.conMinIpsiCong = ga.cueLeft;
contrasts.conMinIpsiCong.individual = ga.conMinIpsiCongruentInvalid.individual- ...
    ga.conMinIpsiCongruentValid.individual;

contrasts.conMinIpsiIncong = ga.cueLeft;
contrasts.conMinIpsiIncong.individual = ga.conMinIpsiIncongruentInvalid.individual - ...
    ga.conMinIpsiIncongruentValid.individual

%contrast: con-ipsi valid / invalid

contrasts.conMinIpsi_V = ga.cueLeft;
contrasts.conMinIpsi_V.individual = ga.conMinIpsiIncongruentValid.individual - ...
    ga.conMinIpsiCongruentValid.individual;


%%
%plotting data figures
%VALID - congruency effect - timelocked to target
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = 'Cz'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 20]
figure; ft_singleplotER(cfg,ga.incong_v,ga.cong_v,contrasts.incongruency_V) % collapsed across right,left
set(gca,'YDir','reverse');
legend('incongruent', 'congruent','incongruent-congruent','location','best')
xlabel('Time (ms)')
ylabel('uV')
title('Cued')
plot2svg('valid_congruency.svg')

%INVALID congruency effect - timelocked to target 
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = 'Cz'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 20]
figure; ft_singleplotER(cfg,ga.incong_i,ga.cong_i,contrasts.incongruency_I) % collapsed across right,left
set(gca,'YDir','reverse');
legend('incongruent', 'congruent','incongruent-congruent','location','best')
xlabel('Time (ms)')
ylabel('uV')
title('Uncued')
plot2svg('invalid_congruency.svg')

% valid congruency - invalid congruency
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = 'Cz'
cfg.xlim = [-.5 1.5]
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 20]
figure; ft_singleplotER(cfg,contrasts.incongruency_I,contrasts.incongruency_V);
xlabel('Time (ms)')
ylabel('uV')
set(gca,'YDir','reverse');
title('Congruency Effect: Difference waves - CZ');
legend('Invalid','Valid','location','best')
plot2svg('inval_cong-vali_cong.svg')



cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.baseline = [-0.2 0]
cfg.xlim = [-0.5 1.5]
%cfg.ylim = [-10 10]
figure; ft_multiplotER(cfg,ga.cueLeft,ga.cueRight);
legend('cue left', 'cue right','location','northwest')
xlabel('Time (ms)')
ylabel('uV')
title('Attentional Preparation')


% Plot the CNV
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.baseline = [-0.2 0]
cfg.channel = 'Cz';
cfg.xlim = [-0.5 1.5]
%cfg.ylim = [-10 10]
figure; ft_singleplotER(cfg,ga.cueLeft,ga.cueRight);
set(gca,'YDir','reverse');
legend('cue left', 'cue right','location','northwest')
xlabel('Time (ms)')
ylabel('uV')
title('Attentional Preparation')
%plot2svg('cnv.svg')

%%
%plot figures


%congruency
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 20]
set(gca,'YDir','reverse');
figure; ft_multiplotER(cfg,ga.incongruent,ga.congruent,contrasts.incongruency)
legend('incongruent','congruent','location','best');
title('congruency effect')
plot2svg('congruency_multi.svg')

cfg = [];
%cfg.layout = layoutmw64;
cfg.channel = 'Cz'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_singleplotER(cfg,ga.incongruent,ga.congruent,contrasts.incongruency)
set(gca,'YDir','reverse');
legend('incongruent','congruent','incongruent-congruent','location','best');
title('congruency - Cz');
plot2svg('congruency.svg')



%plot singel channel 
%VALID - INVALID
%time locked to target
cfg = [];
%cfg.layout = layoutmw64;
cfg.channel = 'Cz'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_singleplotER(cfg,ga.invalid,ga.valid,contrasts.validity)
set(gca,'YDir','reverse');
legend('Uncued,','Cued','Uncued-Cued','location','best');
title('Cue Condition - Cz');
plot2svg('validity.svg')


%time locked to target
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels    = 'yes'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_multiplotER(cfg,ga.invalid,ga.valid,contrasts.validity)
title('Validity');
plot2svg('validity_multi.svg')



%%
%incongruent invalid - incongruent valid
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels    = 'yes'
cfg.baseline = [-0.2 0]
cfg.ylim = [-5 15]
cfg.xlim = [-0.2 1]
figure; ft_multiplotER(cfg,ga.incong_i,ga.incong_v,contrasts.validity_I)
title('Cue effect for incongruent trials');
plot2svg('validity-incong_multi.svg')

cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = '36'
cfg.xlim = [-.5 1]
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_singleplotER(cfg,ga.incong_i,ga.incong_v,contrasts.validity_I)
set(gca,'YDir','reverse');
title('Cue effect for Incongruent Trials - Cz')
xlabel('Time (ms)')
ylabel('uV')
legend('Uncued','Cued','Uncued-Cued','location','best')
plot2svg('validity_incong.svg')


%congruent invalid - congruent valid
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels    = 'yes'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_multiplotER(cfg,ga.cong_i,ga.cong_v,contrasts.validity_C)
title('Cue effect for Congruent trials');
plot2svg('validity_cong_multi.svg')

cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = '36'
cfg.xlim = [-.5 1]
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_singleplotER(cfg,ga.cong_i,ga.cong_v,contrasts.validity_C)
set(gca,'YDir','reverse');
title('Cue effect for Congruent Trials - Cz')
xlabel('Time (ms)')
ylabel('uV')
legend('Uncued','Cued','Uncued-Cued','location','best')
plot2svg('validity_cong.svg')








%%
%plots

%incongruency effect
load  layoutmw64.mat;
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels    = 'yes'
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 10]
figure; ft_multiplotER(cfg,ga.incongruent,ga.congruent,contrasts.incongruency)


%plot singel channel for congruency - 
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = 'Cz'
cfg.xlim = [-.5 1]
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 15]
figure; ft_singleplotER(cfg,ga.incongruent,ga.congruent, contrasts.incongruency)
set(gca,'YDir','reverse');
title('Incongruency Effect - Cz')
xlabel('Time (ms)')
ylabel('uV')
legend('incongruent','congruent','incongruent - congruent','location','best')
%plot2svg('congruency.svg')


% DIFFERENCE WAVES valid congruency - invalid congruency
load  layoutmw64.mat
cfg = [];
cfg.layout = layoutmw64;
cfg.showlabels  = 'yes'
cfg.channel = 'Cz'
cfg.xlim = [-.5 1]
cfg.baseline = [-0.2 0]
cfg.ylim = [-10 10]
figure; ft_singleplotER(cfg,contrasts.incongruency_I,contrasts.incongruency_V);
%could add a difference wave for these 
xlabel('Time (ms)')
ylabel('uV')
title('incongruent-congruent for valid, invalid - Cz');
set(gca,'YDir','reverse');
legend('invalid','valid','location','best')

