function data=ft_remove_ems_v2(data,booleans)
% Simple micro to remove the mean

% 1) create a nice 3d vector
for ii=1:length(data.trial)
    ind(ii,:,:)=data.trial{ii};
end
% 2) now for each channel and contrast remove that mean
fn=fieldnames(booleans); c=double(zeros(1,length(booleans.(fn{1}))));
for icond=1:length(fn)
    ind_avg=mean(ind(booleans.(fn{icond}),:,:));
    c=c+booleans.(fn{icond});
    for ii=find(booleans.(fn{icond}))
        ind(ii,:,:)=ind(ii,:,:)-ind_avg;
        % 3) now reshape this back into data.trial
        data.trial{ii}=squeeze(ind(ii,:,:));
    end
end

if max(c)>2
    display('In ft_remove_ems_v2.m');
    display('Your chosen conditions are not independent and thus things broke');
    keyboard;
end

