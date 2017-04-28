function granger_calc(tt)
global g;

if g.Model_Opt==1
    % Compute best model order
    [AIC,BIC,moAIC,moBIC]=tsdata_to_infocrit(g.X,g.momax,g.icregmode);
    g.BIC(tt,:)=BIC;
    g.AIC(tt,:)=AIC;
    % 
    % figure(1); clf;
    % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    % title('Model order estimation');    
    if     strcmpi(g.str_morder,'actual')
        g.smorder(tt) = amo;
        fprintf('\nusing actual model order = %d\n',g.morder);
    elseif strcmpi(g.str_morder,'AIC')
        g.smorder(tt) = moAIC;
        fprintf('\nusing AIC best model order = %d\n',g.morder);
    elseif strcmpi(g.str_morder,'BIC')
        g.smorder(tt) = moBIC;
        fprintf('\nusing BIC best model order = %d\n',g.morder);
    else
        fprintf('\nusing specified model order = %d\n',g.morder);
    end
    return;
end
%-------------------------%
% Vectorized Autoregression
%-------------------------%
[A,SIG]=tsdata_to_var(g.X,g.morder,'LWR');
assert(~isbad(A),'VAR estimation failed');
% A -> regression coeffs
% SIG -> residuals covariance matrix (key-output of many models)
%-------------------------%
% Checking...
%-------------------------%
[G,info] = var_to_autocov(A,SIG);
var_info(info,true); % report results (and bail out on error)

% Temporal
g.time.F = autocov_to_pwcgc(G);
assert(~isbad(g.time.F,false),'GC calculation failed');
    
% g.time.pval = mvgc_pval(g.time.F,g.morder,g.nobs,g.ntrials,1,1,g.nvars-2,g.tstat); % take careful note of arguments!
% sig  = significance(g.time.pval,g.alpha,g.mhtc);

% figure(2); clf;
% subplot(1,3,1);
% plot_pw(g.time.F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(g.time.pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])

% cd = mean(F(~isnan(F)));
% fprintf('\ncausal density = %f\n',cd);
if g.freq_calc==1
    
% Spectral
g.freq.f = autocov_to_spwcgc(G,g.fres);
assert(~isbad(g.freq.f,false),'spectral GC calculation failed');
% figure(3); clf;
% plot_spw(g.freq.f,g.fs);

% g.freq.pval = mvgc_pval(g.freq.f,g.morder,g.nobs,g.ntrials,1,1,g.nvars-2,g.tstat);
% plot_spw(pval,fs);

% fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
% Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
% mad = maxabs(F-Fint);
% madthreshold = 1e-5;
% if mad < madthreshold
%     fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
% else
%     fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
% end
end
