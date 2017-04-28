function pull_peaks_bar(X1,f1,f2,s,Y)

Nsubj=length(unique(X1(:,s)));
for ii=1:Nsubj
    subjI=X1(:,s)==ii;
    for jj=1:f1.s
        if f1.d==1
            % Low to high
            A=unique(X1(:,f1.I));
            A(A==0)=[];
            f1I=X1(:,f1.I)==A(jj);
        else
            A=mean(X1(:,f1.I));
            if jj==1, f1I=X1(:,f1.I)<=A;
            elseif jj==2, f1I=X1(:,f1.I)>A;
            end
        end
        clear A;
        for kk=1:f2.s
            if f2.d==1
                % Low to high
                A=unique(X1(:,f2.I));
                A(A==0)=[];
                f2I=X1(:,f2.I)==A(jj);
            else
                A=mean(X1(:,f2.I));
                if kk==1, f2I=X1(:,f2.I)<=A;
                elseif kk==2, f2I=X1(:,f2.I)>A;
                end
            end
            clear A;
            ca(ii,jj,kk)=mean(Y([subjI & f1I & f2I]));
        end
        
    end
end

figure(1);
bar(squeeze(mean(ca)))
set(gca,'XTickLabel',f2.n);
legend(f1.n);

