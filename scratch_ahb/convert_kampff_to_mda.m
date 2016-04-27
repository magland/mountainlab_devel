% script to convert Kamppf juxta paired MEA data to MDA, and extract
% the JC firing events.
% Barnett 4/27/16
% See READ ME.docx and SOFTWARE.docx (python file-reading code) in each subdir.
clear

dir = '2014_11_25_Pair_3_0'; %%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen([dir,'/amplifier2014-11-25T23_00_08.bin'],'r');       % raw MEA
Y = fread(fid,'uint16');  % size not given
fclose(fid);
M=32;   % # ch
N = numel(Y)/M
Y = reshape(Y,[M N]);
Y = (Y-32768)*0.195;    % uV
%figure; plot(Y(1,1:1e4))
writemda32(Y,[dir,'/raw.mda']);  % is 2 GB

fid=fopen([dir,'/adc2014-11-25T23_00_08.bin']);        % ADC = juxta
J = fread(fid,'uint16');
fclose(fid);
MJ = 8;   % # adc ch
N = numel(J)/MJ
J = reshape(J,[MJ N]);
used_channel = 0; J = J(used_channel+1,:);  % to 1-indexed channel #
J = J * (10/65536/100) * 1e6;     % uV
%spikespy(J);
writemda32(J,[dir,'/juxta.mda']);
mJ = mean(J);
times = find(diff(J>mJ+(max(J)-mJ)/2)==1);  % trigger on half-way-up-going
labels = 1+0*times;
writemda64([0*times;times;labels],[dir,'/truefirings.mda']);

% elec coords (x,y)
ord = [31 24 7 1 21 10 30 25 6 15 20 11 16 26 5 14 19 12 17 27 4 8 18 13 23 28 3 9 29 2 22]; % ordering across, pasted from PDF file Map_32electrodes.pdf, apart from the top 0.
ord = ord+1;   % 1-indexed
x=nan(32,1); y=x;
x(1) = 0; x(ord(1:3:end))=0;
x(ord(2:3:end))=-sqrt(3)/2; x(ord(3:3:end))=+sqrt(3)/2;
y(1) = 0; y(ord(1:3:end))=-1:-1:-11;
y(ord(2:3:end))=-1.5:-1:-10.5; y(ord(3:3:end))=-1.5:-1:-10.5;
figure; plot(x,y,'k.'); hold on; title('0-indexed electrode locations');
for m=1:M, text(x(m),y(m),sprintf('%d',m-1)); end, axis equal

a = zeros(M,M); % generic build adjacency matrix via dist
for m=1:M, a(m,find(sqrt((x-x(m)).^2+(y-y(m)).^2)<1.5)) = 1; end
% check it
for m=1:M, for k=1:M, if a(m,k), plot([x(m) x(k)],[y(m) y(k)],'r-'); end, end, end
writemda32(a,[dir,'/elecadjacencymatrix.mda']);


dir = '2014_03_26_Pair_2_0'; %%%%%%%%%%%%%%%%%%%%%%%%% another dataset

fid=fopen([dir,'/amplifier2014-03-26T05_11_53.bin'],'r');       % raw MEA
Y = fread(fid,'float32');  % size not given
fclose(fid);
M=32;   % # ch
N = numel(Y)/M
Y = reshape(Y,[M N]);   % already in uV (no subtraction or scaling needed)
writemda32(Y,[dir,'/raw.mda']);  % is 2 GB

fid=fopen([dir,'/adc2014-03-26T05_11_53.bin']);        % ADC = juxta
J = fread(fid,'uint32');
fclose(fid);
MJ = 8;   % # adc ch
N = numel(J)/MJ
J = reshape(J,[MJ N]);
used_channel = 1; J = J(used_channel+1,:);
J = J * (10/65536/1000) * 1e6;     % uV
writemda32(J,[dir,'/juxta.mda']);              % not really needed
mJ = mean(J);
times = find(diff(J>mJ+(max(J)-mJ)/2)==1);  % trigger on half-way-up-going
labels = 1+0*times;
%spikespy({J,times,labels});
writemda64([0*times;times;labels],[dir,'/truefirings.mda']);

writemda32(a,[dir,'/elecadjacencymatrix.mda']);  % same MEA layout

% note: 30kHz sample rate
