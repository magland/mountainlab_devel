function test_noise_cutoff

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/../franklab/experiment1/output',mfile_path);

clips=readmda(sprintf('%s/clips.mda',path0));
cii=readmda(sprintf('%s/clips_index.mda',path0));
K=length(cii);
[M,T,NC]=size(clips);
cii=[cii,NC];

spikes=22;

for k=spikes
    ii=(cii(k)+1):(cii(k+1));
    clips0=clips(:,:,ii);
    template0=mean(clips0,3);
    maxes0=max(abs(template0),[],2);
    [~,load_channel]=max(maxes0);
    tmp=squeeze(max(abs(clips0(load_channel,:,:)),[],2));
    figure; hist(tmp,1000);
end;

end
