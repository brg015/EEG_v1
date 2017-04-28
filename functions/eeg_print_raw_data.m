function save_data=eeg_print_raw_data(iSubj)
global RUN;

data=excel_reader(fullfile(RUN.dir.beh,['subject' RUN.dir.subjects{iSubj} '.csv']));
for ii=1:length(data), head{ii}=data{ii}.header{1}; end
% Grab some relevent fields
catch_T=find(strcmp(head,'Catch'));
E1_catch_R=find(strcmp(head,'E1_catch_R'));
E1_catch_hit=find(strcmp(head,'E1_catch_hit'));
R1_hc=find(strcmp(head,'R1_hc'));
R1_hit=find(strcmp(head,'R1_hit'));
R1_miss=find(strcmp(head,'R1_miss'));
R1_cr=find(strcmp(head,'R1_cr'));
R1_fa=find(strcmp(head,'R1_fa'));
R2_hc=find(strcmp(head,'R2_hc'));
R2_hit=find(strcmp(head,'R2_hit'));
R2_miss=find(strcmp(head,'R2_miss'));
R2_cr=find(strcmp(head,'R2_cr'));
R2_fa=find(strcmp(head,'R2_fa'));
R2_s=find(strcmp(head,'R2_s'));
R1R2_on=find(strcmp(head,'R1R2_ON'));
R1R2_no=find(strcmp(head,'R1R2_NO'));
R1R2_oo=find(strcmp(head,'R1R2_OO'));
R1R2_nn=find(strcmp(head,'R1R2_NN'));
c1=find(strcmp(head,'ID'));
descrip_idx=[E1_catch_R E1_catch_hit R1_hc R1_hit R1_miss R1_cr R1_fa ...
    R2_hc R2_hit R2_miss R2_cr R2_fa R1R2_on R1R2_no R1R2_oo R1R2_nn];
descrip={'E1_R' 'E1_hit' 'R1_hc' 'R1_hit' 'R1_miss' 'R1_cr' 'R1_fa' ...
    'R2_hc' 'R2_hit' 'R2_miss' 'R2_cr' 'R2_fa' 'R1R2_on' 'R1R2_no' 'R1R2_oo' 'R1R2_nn'};
        
for ii=1:length(data{1}.col);
    for jj=1:length(descrip_idx)
        descrip_mat(ii,jj)=str2double(data{descrip_idx(jj)}.col{ii});
    end
    code_mat(ii)=str2double(data{c1}.col{ii});
end
        
% Now let's save this output
Tcounts=nansum(descrip_mat);
for ii=1:length(Tcounts),
    save_data{ii}.header=descrip{ii};
    save_data{ii}.col(1)=Tcounts(ii);
end
        
[~,type,~]=fileparts(RUN.dir.beh);
write_struct(save_data,fullfile(RUN.dir.QAL,...
    'Subject',RUN.dir.subjects{iSubj},[type '_raw_epoch_count.csv']));