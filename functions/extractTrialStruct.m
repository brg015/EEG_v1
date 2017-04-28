function trialStruct = extractTrialStruct(EEG,iSubj)
%=========================================================================%
% BVD - BRG edits Summer 2014
%=========================================================================%
global RUN; eeglab -redraw;
%-------------------------%
% initialize events
%-------------------------% 
% This recently broke, but it isn't even used anyways, just commenting it
% out for now
% for i= 1:length(EEG.event)
%     events(i) = str2double(EEG.event(i).type);
%     times(i) = EEG.event(i).latency ;
% end

switch RUN.dir.study
    case 'SEFER'
        %-------------------------%
        % Addin behave data
        %-------------------------%
        data=excel_reader(fullfile(RUN.dir.beh,['subject' RUN.dir.subjects{iSubj} '.csv']));
        % Need to add in COCA_lemma stuff to specifc trials
        %=================================================================%
        % Lemma Micro
        %=================================================================%
        COCA=excel_reader(fullfile(RUN.template.COCA));
        C_ID=cell2num(COCA{1}.col);
        ID=cell2num(data{1}.col);
        % Both Sorted by ID so this is kinda easy
        if sum(C_ID-ID)~=0
            fprintf('Oh shit, COCA, don''t match\n');
            fprintf('Go to extractTrailStruct');
            keyboard;
        end
        clear ID C_ID
        Dom=cell2num(COCA{13}.col);
        High_Freq=cell2num(COCA{8}.col);
        Low_Freq=cell2num(COCA{9}.col);
        High_L_Freq=cell2num(COCA{19}.col);
        Low_L_Freq=cell2num(COCA{20}.col);
        Nwords=cell2num(COCA{21}.col);
        CAT=unique(COCA{11}.col);
        SIZE=cell2num(COCA{15}.col);
        for jj=1:length(CAT),
            Catg(strcmp(COCA{11}.col,CAT{jj}))=jj;
        end
        %=================================================================%

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
            'R2_hc' 'R2_hit' 'R2_miss' 'R2_cr' 'R2_fa' 'R1R2_on' 'R1R2_no' 'R1R2_oo' 'R1R2_nn' ...
            'HighFreq' 'LowFreq' 'TrialID' 'High_L_Freq' 'Low_L_Freq' 'Compound' 'Catg' 'Domain' 'Size'};
    
        for ii=1:length(data{1}.col);
            for jj=1:length(descrip_idx)
                descrip_mat(ii,jj)=str2double(data{descrip_idx(jj)}.col{ii});
            end
            code_mat(ii)=str2double(data{c1}.col{ii});
        end
        
        trialStruct.descrip=descrip;
        % Epochs can be pre, stim or R
        c=0;

        for ii=1:length(EEG.epoch)
            trialStruct.trial(ii)=c;
            trailStruct.epoch(ii)=EEG.epoch(ii);   
            trialStruct.phase(ii)=0;
            % Find the event code with latency of 0
            try
                eventtype_idx=find(cell2mat(EEG.epoch(ii).eventlatency)==0);
                if length(eventtype_idx)>1
                    display(['WARNING: overlapping indices in epoch ' num2str(ii)]);
                    display('...Something may be wrong');
                    eventtype=str2double(EEG.epoch(ii).eventtype(...
                        find(~strcmp(EEG.epoch(ii).eventtype(eventtype_idx),'999'))));
                else
                    eventtype=str2double(EEG.epoch(ii).eventtype(eventtype_idx));
                end
            catch err
                keyboard
            end
            
            if (eventtype > 10000 && eventtype<40000)  
                % We never get into here anymore, as we're looking at stim
                % response only
                offset_v=double(int8(eventtype/10000))*10000;
                c=c+1;
                rProfile=[];
                trialStruct.event(ii)=1; % pre
                ID(c)=find(code_mat==(eventtype-offset_v)); 
                rProfile=descrip_mat(ID(c),:); 
                descrip_rep(c,:)=rProfile; 
                trialStruct.trial(ii)=c;   
            elseif (eventtype > 40000)
                % Offset_v accounts for the addition of 30000 to the code
                offset_v=double(int8((eventtype-30000)/10000))*10000;
                c=c+1;
                rProfile=[];
                trialStruct.event(ii)=2; % stim present
                % Now we need the rProfile as well and ID codes
                ID(c)=find(code_mat==(eventtype-offset_v-30000)); 
                rProfile=descrip_mat(ID(c),:); 
                descrip_rep(c,1:16)=rProfile; 
                trialStruct.trial(ii)=c;   
                % Add in freq info
                TRIAL_id=str2num(data{1}.col{ID(c)});
                rProfile(17)=High_Freq(ID(c));
                rProfile(18)=Low_Freq(ID(c));
                rProfile(19)=TRIAL_id;
                rProfile(20)=High_L_Freq(ID(c));
                rProfile(21)=Low_L_Freq(ID(c));
                rProfile(22)=Nwords(ID(c));
                rProfile(23)=Catg(ID(c));
                rProfile(24)=Dom(ID(c));
                rProfile(25)=SIZE(ID(c));
                descrip_rep(c,17)=High_Freq(ID(c)); % High Freq
                descrip_rep(c,18)=Low_Freq(ID(c));  % Low Freq
                descrip_rep(c,19)=TRIAL_id;
                descrip_rep(c,20)=High_L_Freq(ID(c));
                descrip_rep(c,21)=Low_L_Freq(ID(c));   
                descrip_rep(c,22)=Nwords(ID(c)); 
                descrip_rep(c,23)=Catg(ID(c));
                descrip_rep(c,24)=Dom(ID(c));
                descrip_rep(c,25)=SIZE(ID(c));
            else
                trialStruct.event(ii)=3; % response
            end
            trialStruct.rProfile(ii,:)=rProfile;
            trialStruct.phase(ii)=offset_v/10000;
        end  
        trialStruct.eProfile=descrip_rep;
        %=================================================================%
        % Now alter trials to remove catch trials from hits or misses
        % But preserve info on trial type*
        %=================================================================%
        % Probably just wanna kill eProfile at some point, but for now,
        % let's save it and move on.
        
        % Catch FIX - Changes behave ONLY
        bI=trialStruct.rProfile(:,1)==1;
        trialStruct.rProfile(bI,3:16)=0;
        clear bI;
        bI=trialStruct.eProfile(:,1)==1;
        trialStruct.eProfile(bI,3:16)=0;
        clear bI;
        
        % Compound FIX - Changes everything
        bI=trialStruct.rProfile(:,22)>1;
        % Save trial ID, but kill others
        trialStruct.rProfile(bI,3:18)=0;
        trialStruct.rProfile(bI,20:21)=0; 
        trialStruct.rProfile(bI,23:25)=0;
        
        clear bI;
        bI=trialStruct.eProfile(:,22)>1;
        % Save trial ID, but kill others
        trialStruct.eProfile(bI,3:18)=0;
        trialStruct.eProfile(bI,20:21)=0; 
        trialStruct.eProfile(bI,23:25)=0; 
        %=================================================================%
        
        % And set NaNs to zeros
        [ii,jj]=find(isnan(trialStruct.rProfile));
        for aa=1:length(ii), trialStruct.rProfile(ii(aa),jj(aa))=0; end
        clear ii jj;
        
        [ii,jj]=find(isnan(trialStruct.eProfile));
        for aa=1:length(ii), trialStruct.eProfile(ii(aa),jj(aa))=0; end
        
end



