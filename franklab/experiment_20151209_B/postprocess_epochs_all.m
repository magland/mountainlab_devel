function postprocess_epochs_all

mfile_path=fileparts(mfilename('fullpath'));

strings={'s1','r1','s2','r2','s3','r3'};
Ns=[];
for j=1:length(strings)
    str=strings{j};
    disp(str);
    pre2a=readmda([mfile_path,sprintf('/output_%s/pre2.mda',str)]);
    Ns(j)=size(pre2a,2);
end

firings=readmda([mfile_path,'/output_epochs_all/firings.mda']);

times=firings(2,:);
labels=firings(3,:);

for j=1:length(strings)
    str=strings{j};
    disp(str);
    indsA=find((times>sum(Ns(1:j-1)))&(times<=sum(Ns(1:j))));
    firingsA=firings(:,indsA);
    firingsA(2,:)=firingsA(2,:)-sum(Ns(1:j-1));
    mv.raw=[mfile_path,sprintf('/output_%s/pre2.mda',str)];
    mv.firings=firingsA;
    mv.sampling_freq=30000;
    ms_mountainview(mv);
end

end