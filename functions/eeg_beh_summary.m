function eeg_beh_summary()
global RUN;
% Make a mini summary of included events
get_events={'E_r' 'E_hit' 'R1hit' 'R1miss' 'R2hit' 'R2miss' 'R1R2oo' 'R1R2nn' 'R1R2on' 'R1R2no'};

% Setup group output
for jj=2:length(get_events), 
    g_data{jj}.header=get_events{jj-1}; 
    g_data_N{jj}.header=get_events{jj-1};
end
g_data{1}.header='Subject';
g_data_N{1}.header='Subject';

ttypes={'pre','stm','rsp'};

for iSubj=1:length(RUN.dir.subjects)  
    % All those save files
    [~,save_prefix,~]=fileparts(['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]);
    QA_file=fullfile(RUN.dir.QAL,'Subject',RUN.dir.subjects{iSubj},[save_prefix '_N.csv']);
    QA_file2=fullfile(RUN.dir.QAL,'Subject',RUN.dir.subjects{iSubj},[save_prefix '_N_basic.csv']);
    if exist(QA_file,'file'), data=excel_reader(QA_file); else continue; end
    
    % Setup Out Data to save
    out_data{1}.header='Counts';
    for jj=2:length(get_events)+1; 
        out_data{jj}.header=get_events{jj-1}; 
    end
    
    % Collect raw data
    save_data=eeg_print_raw_data(iSubj);
    % Remove pesky _
    for jj=1:length(save_data)
        kI=strfind(save_data{jj}.header,'_');
        SDH{jj}=...
            save_data{jj}.header(setdiff(1:length(save_data{jj}.header),kI));
    end
    SDH{1}='E_r';
    SDH{2}='E_hit';
    
    % Add info to out_data
    out_data{1}.col{1}='Raw';
    for jj=1:length(get_events)
        X=strcmp(get_events{jj},SDH);
        V=save_data{X}.col;
        out_data{jj+1}.col(1)=V;
    end
    
    out_data{1}.col{2}='pre';
    out_data{1}.col{3}='stm';
    out_data{1}.col{4}='rsp';
    for jj=1:length(ttypes)
        I=strmatch(ttypes{jj},data{1}.col);
        Tmatch=data{1}.col(I);
        Nmatch=cell2num(data{2}.col(I));
        for kk=1:length(get_events)
            X=strsearch(get_events{kk},Tmatch);
            V=Nmatch(X);
            out_data{kk+1}.col(jj+1)=V;
        end
    end
    write_struct(out_data,QA_file2);
    %==============%
    % Group Output:
    %==============%
    g_data{1}.col{iSubj}=RUN.dir.subjects{iSubj};
    g_data_N{1}.col{iSubj}=RUN.dir.subjects{iSubj};
    for jj=2:length(get_events),
        % Ratio of stim to include events
        g_data_N{jj}.col{iSubj}=out_data{jj}.col(3);
        g_data{jj}.col(iSubj)=out_data{jj}.col(3)/out_data{jj}.col(1);
    end
end

[~,hc_set,~]=fileparts(RUN.dir.beh);
write_struct(g_data_N,fullfile(RUN.dir.QAL,'Group',[hc_set '_N.csv']));
write_struct(g_data,fullfile(RUN.dir.QAL,'Group',[hc_set '_per.csv']));