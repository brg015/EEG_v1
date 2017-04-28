function granger_plot(sv2to1,sv1to2,sFv1to2,sFv2to1,ti,save_dir,save_file,bb)
global g;
[~,nam,~]=fileparts(save_file);
lam=sfreqs(125,g.fs); % Hardcoded
sf1=fullfile(save_dir,[nam '_Timeseries.jpg']);
sf2=fullfile(save_dir,[nam '_BetaBand.jpg']);
sf3=fullfile(save_dir,[nam '_SpectMaps.jpg']);
sf4=fullfile(save_dir,[nam '_SpectT1.jpg']);
sf5=fullfile(save_dir,[nam '_SpectT2.jpg']);
    
if bb==1
    figure(1);
    subplot(2,2,1);
    plot(ti,mean(sv1to2{1}),'b','linewidth',3); hold on;
    plot(ti,mean(sv1to2{2}),'r','linewidth',3); 
    legend({'Condition 1','Condition 2'});
    title('X to Y'); grid;
    subplot(2,2,2);
    plot(ti,mean(sv2to1{1}),'b','linewidth',3); hold on;
    plot(ti,mean(sv2to1{2}),'r','linewidth',3); 
    legend({'Condition 1','Condition 2'});
    title('Y to X'); grid;
    subplot(2,2,3);
    [~,~,~,s1]=ttest(sv1to2{1}-sv1to2{2});
    plot(ti,s1.tstat,'k','linewidth',3); 
    title('X to Y (ttest)'); grid;
    subplot(2,2,4);
    [~,~,~,s2]=ttest(sv2to1{1}-sv2to1{2});
    plot(ti,s2.tstat,'k','linewidth',3); 
    title('Y to X (ttest)'); grid;
    set(gcf,'position',[0 0 1280 1024]); 
    set(gcf,'color','w');
    export_fig(sf1);
end

if bb==2
    
    FI=(lam>=16 & lam<=18);
    for ii=1:2
        A{ii}=squeeze(mean(sFv1to2{ii}(:,:,FI),3));
        B{ii}=squeeze(mean(sFv2to1{ii}(:,:,FI),3));
    end
    figure(2);
    subplot(2,2,1);
    plot(ti,mean(A{1}),'b','linewidth',3); hold on;
    plot(ti,mean(A{2}),'r','linewidth',3); 
    legend({'Condition 1','Condition 2'});
    title('X to Y'); grid;
    subplot(2,2,2);
    plot(ti,mean(B{1}),'b','linewidth',3); hold on;
    plot(ti,mean(B{2}),'r','linewidth',3); 
    legend({'Condition 1','Condition 2'});
    title('Y to X'); grid;
    subplot(2,2,3);
    [h1,p1,~,s1]=ttest(A{1}-A{2});
    plot(ti,s1.tstat,'k','linewidth',3); 
    title('X to Y (ttest)'); grid;
    subplot(2,2,4);
    [p2,~,~,s2]=ttest(B{1}-B{2});
    plot(ti,s2.tstat,'k','linewidth',3); 
    title('Y to X (ttest)'); grid;
    set(gcf,'position',[0 0 1280 1024]); 
    set(gcf,'color','w');
    export_fig(sf2);

    try
    FI=(lam>=1 & lam<=60);
    for ii=1:2
        A{ii}=sFv1to2{ii}(:,:,FI);
        B{ii}=sFv2to1{ii}(:,:,FI);
    end

    figure(3);
    subplot(2,3,1);
    imagesc(squeeze(mean(A{1}))); title('X to Y Condition 1');
    subplot(2,3,2);
    imagesc(squeeze(mean(A{2}))); title('X to Y Condition 2');
    subplot(2,3,3);
    imagesc(squeeze(mean(A{1}))-squeeze(mean(A{2}))); title('X to Y Condition 1-2');
    subplot(2,3,4);
    imagesc(squeeze(mean(B{1}))); title('Y to X Condition 1');
    subplot(2,3,5);
    imagesc(squeeze(mean(B{2}))); title('Y to X Condition 2');
    subplot(2,3,6);
    imagesc(squeeze(mean(B{1}))-squeeze(mean(B{2}))); title('Y to X Condition 1-2');
    set(gcf,'position',[0 0 1280 1024]); 
    set(gcf,'color','w');
    export_fig(sf3);

    figure(4);
    [~,~,~,s1]=ttest(A{1}-A{2});
    imagesc(squeeze(s1.tstat)); title('X to Y'); colorbar;
    set(gcf,'position',[0 0 1280 1024]); 
    set(gcf,'color','w');
    export_fig(sf4);

    figure(5);
    [~,~,~,s2]=ttest(B{1}-B{2});
    imagesc(squeeze(s2.tstat)); title('Y to X'); colorbar;
    set(gcf,'position',[0 0 1280 1024]); 
    set(gcf,'color','w');
    export_fig(sf5);
    catch err
        display('Bad Freq. Info');
    end

end

close all;
