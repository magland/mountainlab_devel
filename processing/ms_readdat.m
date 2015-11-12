function X=ms_readdat(path,opts)

channels=opts.channels;
num_channels=opts.num_channels;
timepoints=opts.timepoints;
dtype=opts.dtype;

if ~isempty(timepoints)
    min_timepoint=min(timepoints);
    max_timepoint=max(timepoints);
else
    min_timepoint=1;
end;

F=fopen(path,'rb');
if (min_timepoint>1)
    fread(F,[num_channels,min_timepoint-1],dtype);
end;
if (~isempty(timepoints))
    X=fread(F,[num_channels,max_timepoint-min_timepoint+1],dtype);
else
    X=fread(F,[num_channels,inf],dtype);
end;
fclose(F);

if (~isempty(timepoints))
    X=X(channels,timepoints-min_timepoint+1);
else
    X=X(channels,:);
end;

end

