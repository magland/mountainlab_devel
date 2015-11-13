function CC_stats

CC=mda_to_cross_correlograms(readmda('example1_output/cross-correlograms.mda'));
%view_cross_correlograms(CC,0);

n1=50;
n2=100;

show_CC_stats(

K=size(CC,1);
AA=zeros(K,K);
for k1=1:K
    for k2=1:K
        tmp=CC{k1,k2};
        kk=length(find(abs(tmp)<n1));
        nn=length(find(abs(tmp)<n2));
        if (nn>0)
            pp=n1/n2;
            stdev0=sqrt(nn*pp*(1-pp));
            num_sigmas=(kk-nn*pp)/stdev0;
            if (abs(num_sigmas)>3)
                AA(k1,k2)=log((kk/nn)/pp);
            end;
        end;
        
    end;
end;

figure; imagesc(AA);

end