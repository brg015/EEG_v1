function classify_behave(behav)

switch behav.hc
    case 0, post_fix='';
    case 1, post_fix='_hc';
end

group_file=fullfile(behav.save_dir,['pre' post_fix],'subject00.csv');
template_file=fullfile('D:\Data','SEFER','EEG','template','image_info_v1.csv');

g_data=excel_reader(group_file);          g_head=gethead(g_data);
t_data=excel_reader(template_file);       t_head=gethead(t_data);
for ii=19:37
    g_data{ii}.col=cell2num(g_data{ii}.col);
end
%------------------------
% Sort Data by ID
%------------------------
ID_g=strcmp('ID',g_head);
ID_t=strcmp('ID',t_head);
ID_ga=g_data{ID_g}.col;
ID_ta=t_data{ID_t}.col;
for jj=1:430
    idx=strcmp(ID_ga{jj},ID_ta);
    for ii=1:length(t_data),
        if ii>4
            try
                ID_tas{ii}.col(jj)=str2num(t_data{ii}.col{idx});
            catch err
                ID_tas{ii}.col(jj)=NaN;
            end
        else 
            ID_tas{ii}.col{jj}=t_data{ii}.col{idx};
        end
        if jj==1, ID_tas{ii}.header=t_head{ii}; end
    end
end
clear t_data; t_data=ID_tas; clear ID_tas;
%------------------------
% Setup Output
%------------------------
% Set up filter arrays
cat_set=['all' unique(t_data{4}.col), unique(t_data{3}.col)];
cat_set=cat_set(~cellfun('isempty',cat_set));
cmp_array(1,:)=ones(1,430);
for jj=2:15
    if jj<14, cmp=4; else cmp=3; end
    cmp_array(jj,:)=strcmp(cat_set{jj},t_data{cmp}.col(1:430));
end

% Collect g_data stats
g_array_head={'R1hit' 'R1miss' 'R1cr' 'R1fa' 'R2hit' 'R2miss' 'R2cr' 'R2fa' ...
    'R1R2ON' 'R1R2NO' 'R1R2OO' 'R1R2NN','R1N','R2N','R1R2N'};
c=1;
for ii=[19:22 25:28 33:34 36:37]
    g_array(c,:)=g_data{ii}.col; c=c+1;
end
g_array(13,:)=sum(g_array(1:4,:));
g_array(14,:)=sum(g_array(5:8,:));
g_array(15,:)=sum(g_array(9:12,:));
%------------------------
% Create Output
%------------------------
t_data{19}.header='Corr';
for ii=1:length(cat_set)
   outputfile=fullfile(behav.save_dir,['compiled_measures_' cat_set{ii} post_fix '.csv']); 
   f_array=g_array(:,logical(cmp_array(ii,:)));
   data{1}.header='Trials';
   data{1}.col=g_array_head;
   data{2}.header='Percentage';
   data{2}.col(1:4)=f_array(1:4,:)/f_array(13,:);
   data{2}.col(5:8)=f_array(5:8,:)/f_array(14,:);
   data{2}.col(9:12)=f_array(9:12,:)/f_array(15,:);
   for jj=5:22
       data{jj-2}.header=t_data{jj}.header;
       msrs=t_data{jj}.col(logical(cmp_array(ii,:)));
       
       ck=13;
       good_trials=~or(isnan(msrs),f_array(ck,:)==0);
       for kk=1:4,
           data{jj-2}.col(kk)=f_array(kk,good_trials).*(msrs(good_trials))/f_array(kk,good_trials);
       end
       clear good_trials;
       
       ck=14;
       good_trials=~or(isnan(msrs),f_array(ck,:)==0);
       for kk=5:8
           data{jj-2}.col(kk)=f_array(kk,good_trials).*(msrs(good_trials))/f_array(kk,good_trials);
       end
       clear good_trials;
       
       ck=15;
       good_trials=~or(isnan(msrs),f_array(ck,:)==0);
       for kk=9:12,
           data{jj-2}.col(kk)=f_array(kk,good_trials).*(msrs(good_trials))/f_array(kk,good_trials);
       end
       clear good_trials;
   end
   write_struct(data,outputfile);
   clear outputfile f_array data
end








