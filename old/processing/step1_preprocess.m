function data=step1_preprocess(opts)

channels=opts.channels;

timerA=tic;

fname_in=opts.raw_dat;
fname_out=opts.raw_mda;
channels=opts.channels;

opts_fname=[opts.raw_mda,'.mountainsort.opts'];
if (exist(fname_out,'file'))&&(options_match(opts,opts_fname))
    fprintf('File %s already exists and options match.\n',fname_out);
%     if (~isempty(whos('global','global_X')))
%         global global_X;
%         data.X=global_X;
%     else
        fprintf('Reading %s... ',fname_out);
        if ~isempty(opts.timepoints)
            data.X=readmda_data_beginning(fname_out,max(opts.timepoints));
            data.X=data.X(:,opts.timepoints);
        else
            data.X=readmda(fname_out);
        end;
%        global global_X;
%        global_X=data.X;
        fprintf('\nElapsed: %g seconds',toc(timerA));
        fprintf('\n');
%    end;
    fprintf('%d chan, %g timepoints\n',size(data.X,1),size(data.X,2));
    return;
end;

in_file_ok=1;
if (~exist(fname_in,'file'))
    in_file_ok=0;
else
    dd=dir(fname_in);
%     if (dd.bytes~=7800266736)
%         fprintf('Wrong number of bytes: %d\n',dd.bytes);
%         in_file_ok=0;
%     end;
end;

if (~in_file_ok)
    error('Raw file does not exist: %s.',fname_in);
end;

fprintf('Creating %s...',fname_out');

N=inf;

fprintf('\nReading %s... ',fname_in);
F=fopen(fname_in,'rb');
X=fread(F,[72,N],'int16');
fclose(F);

%Extract group 1
fprintf('Extracting group 1... ');
X=X(channels,:);

fprintf('%d chan, %g timepoints\n',size(X,1),size(X,2));

fprintf('Filtering... ');
X=ss_freqfilter(X,30000,300,10000);

fprintf('Normalizing... ');
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
    X(j,:)=min(20,max(-20,X(j,:)));
end;
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
end;

fprintf('Writing %s... ',fname_out);
writemda(X,fname_out);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

write_options(opts,opts_fname);

data.X=X;

end

function opts_txt=get_options_text(opts)
opts_txt='';
opts_txt=[opts_txt,sprintf('channels=%d\n',opts.channels)];
opts_txt=[opts_txt,sprintf('prewhiten=%d\n',opts.prewhiten)];

opts_txt=[opts_txt,sprintf('timepoints_code=%d,%d,%d,%f\n',...
    min(opts.timepoints),max(opts.timepoints),length(opts.timepoints),mean(opts.timepoints)...
)];
end

function ret=options_match(opts,fname)
opts_txt=get_options_text(opts);
if (~exist(fname,'file')) ret=true; return; end;
file_txt=read_text_file(fname);
ret=strcmp(opts_txt,file_txt);
end

function write_options(opts,fname)
opts_txt=get_options_text(opts);
write_text_file(opts_txt,fname);
end

function txt=read_text_file(fname)
txt=fileread(fname);
end

function write_text_file(txt,fname)
F=fopen(fname,'w');
if (F==-1) 
    error('Unable to open file %s',fname);
    return; 
end;
fprintf(F,'%s',txt);
fclose(F);
end

