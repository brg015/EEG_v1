function [Fv1to2,Fv2to1,v1to2,v2to1]=granger_shell(data,I,ti,td)
% data    => from fieldtrip
% I       => index of trials to keep
% ti      => time interval
% td      => spacing (delta)
global g;
Fv1to2=[]; Fv2to1=[];
v1to2=[];  v2to1=[];

tdata=data;
tdata.trial=tdata.trial(I);
tdata.time=tdata.time(I);
                    
for tt=1:length(ti)
    
    sdisp(['Time Center = ' num2str(ti(tt))],2);

    g.tv=(data.time{1}>=ti(tt)-td/2 & data.time{2}<=ti(tt)+td/2);
    
    granger_preprocess(tdata);
   
    granger_calc(tt);
 
    if g.Model_Opt==0 
        if g.freq_calc==1
            for kk=1:size(g.freq.f,3)
                % Feedback
                Fv1to2(tt,kk)=g.freq.f(2,1,kk);
                % Feedforward
                Fv2to1(tt,kk)=g.freq.f(1,2,kk);
            end
        end
        v1to2(tt)=g.time.F(2,1);
        v2to1(tt)=g.time.F(1,2);
    end

end % Temporal loop
