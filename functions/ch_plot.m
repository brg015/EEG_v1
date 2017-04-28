function ch_plot(ga,ch1,c1,c2)

figure; cs={'b' 'r' 'b' 'k' 'c'};
for ii=1:size(ga.(c1).individual,1)
%     plot(squeeze(ga.(c1).individual(ii,ch1,:)-ga.(c2).individual(ii,ch1,:)),cs{ii}); hold on;
    plot(squeeze(ga.(c1).individual(ii,ch1,:)),cs{ii}); hold on;
end
