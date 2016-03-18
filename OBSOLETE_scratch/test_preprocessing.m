function test_preprocessing

close all;

mfile_path=fileparts(mfilename('fullpath'));
raw_path=sprintf('%s/../franklab/raw/hippocampal/tetrode',mfile_path);
path=sprintf('%s/test_output',mfile_path);
if ~exist(path,'dir') mkdir(path); end;

templates_raw=readmda(sprintf('%s/../franklab/test_hippocampal_01/output/templates_raw.mda',mfile_path));
T0=size(templates_raw,2);
figure;
ms_view_templates(templates_raw);

N0=1e6;
dt=1/30000;
TT=dt*N0;
df=1/TT;

A=zeros(1,N0);
A(1:T0)=templates_raw(2,:,2);
Ahat=fft(A);
figure; 
plot((0:N0-1)*df,abs(Ahat));

%fprintf('Reading raw...\n');
%X=readmda([raw_path,'/dl12_20151208_NNF_r1_tet16_17_subset.mda']);
%size(X)

o_filter.samplefreq=30000;
o_filter.freq_min=300;
o_filter.freq_max=1000;
o_filter.outlier_threshold=500;

if 0
    X=readmda([raw_path,'/dl12_20151208_NNF_r1_tet16_17_subset.mda']);
    %X(2:end,:)=X(2:end,:)-repmat(X(1,:),size(X,1)-1,1);
    %X=X(2:end,:);
    writemda(X,[path,'/raw.mda']);
end;

mscmd_bandpass_filter([path,'/raw.mda'],[path,'/pre1.mda'],o_filter);
mscmd_whiten([path,'/pre1.mda'],[path,'/pre2.mda']);

fprintf('Reading...\n');
pre1=readmda([path,'/pre1.mda']);
pre2=readmda([path,'/pre2.mda']);
pre1=pre1(:,1:N0);
pre2=pre2(:,1:N0);
[M,N]=size(pre1);

for m=1:M
    pre1(m,:)=pre1(m,:)/sqrt(var(pre1(m,:)));
    pre2(m,:)=pre2(m,:)/sqrt(var(pre2(m,:)));
end;

fprintf('spikespy...\n');
spikespy(cat(1,pre1,zeros(1,N),pre2));

X0=pre1(2,:);
X0hat=fft(X0);
figure; plot(df*(0:N-1),abs(X0hat(1:N)));
title('pre1');

X0=pre1(2,:);
X0hat=fft(X0);
figure; plot(df*(0:N-1),abs(X0hat(1:N)));
title('pre2');

end
