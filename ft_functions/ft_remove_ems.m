function data=ft_remove_ems(data,reject)
% Simple micro to remove the mean

% 1) create a nice 3d vector
for ii=1:length(data.trial)
    ind(ii,:,:)=data.trial{ii};
end
% 2) now for each channel remove that mean
ind_avg=mean(ind(~reject,:,:));
for ii=1:length(data.trial)
    ind(ii,:,:)=ind(ii,:,:)-ind_avg;
    % 3) now reshape this back into data.trial
    data.trial{ii}=squeeze(ind(ii,:,:));
end

