function eeg_plot_group_erp(plot3_save,f1,f2,fs,ts,tstep,te,subjects,setS,save_dir,cfg_topo)
Ncol=6;

Nsubj=size(plot3_save.individual,1);
plot3=plot3_save;

c=1; f=1; 
for ii=ts:tstep:te
    cfg_topo.xlim=[ii ii+tstep]; 
    for jj=1:Nsubj+1

        if jj<=Nsubj
            plot3.individual=[];
            plot3.individual=plot3_save.individual(jj,:,:);
        else
            plot3.individual=plot3_save.individual;
        end

        if c==1, 
            save3=fullfile(save_dir,[f1 '-' f2 '_' fs '_group_topo' num2str(f) '_set' num2str(setS) '.png']);
            figure(f); set(gcf,'position',[0 0 1280 1024]); 
            subaxis(Nsubj+2,Ncol,c,jj,'Spacing', 0, 'Padding', 0, 'Margin', 0,'SpacingVert',0.05,'MarginTop',0.05);
            ft_topoplotER(cfg_topo, plot3); 
        end
        subaxis(Nsubj+2,Ncol,c,jj,'Spacing', 0, 'Padding', 0, 'Margin', 0,'SpacingVert',0.05,'MarginTop',0.05);
        ft_topoplotER(cfg_topo, plot3); 
        %============%
        % Plot titles
        %============%
        def_title=[num2str(int32(ii*1000)) ' : ' num2str(int32((ii+tstep)*1000)) 'ms'];
        if c==1
            if jj>Nsubj
                title('Group','FontSize',14,'FontWeight','bold');
            else
                title([subjects{jj}],'FontSize',14,'FontWeight','bold');
            end
        else
            title(def_title,'FontSize',14,'FontWeight','bold');
        end
        axis tight
        axis off 
        %============%
    end
    c=c+1;
    if c==Ncol+1, c=1; 
        export_fig(save3);
        close(f); f=f+1;
    end
end

