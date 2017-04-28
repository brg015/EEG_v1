function base_ga_freq=ga_freq_baseline(ga_freq)

fn=fieldnames(ga_freq);
for ii=1:length(fn)
   if ii==1, 
       freq=ga_freq.(fn{ii}).powspctrm;
   else
       freq=cat(1,freq,ga_freq.(fn{ii}).powspctrm);
   end
   s{ii}=size(
end