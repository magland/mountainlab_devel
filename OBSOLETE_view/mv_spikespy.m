function mv_spikespy(varargin)
%MV_SPIKESPY - Launch mountainview to view a raw/preprocessed dataset with
%times and labels
%
% Syntax:  
% mv_spikespy(X) -- X is the raw timeseries data (2D raw, #channels x #timepoints)
% mv_spikespy({X,T,L}) -- same as previous, but shows label markers
%                length(t)=length(L) = #events
% mv_spikespy(...,struct('sampling_freq',30000)) -- set the sampling frequency
%                in Hz
%
% Other m-files required: ms_mountainview
%
% See also: ms_mountainview

% Author: Jeremy Magland
% Mar 2016; Last revision: 9-Mar-2016

if nargin<1, test_mv_spikespy; end;

opts=struct;
for j=1:length(varargin)
    arg=varargin{j};
    if (isstruct(arg))
        opts=arg;
    end;
end;
%default options
if (~isfield(opts,'sampling_freq')) opts.sampling_freq=0; end;
if (isfield(opts,'sampfreq')) opts.sampling_freq=opts.sampfreq; end;

raw_arrays={};
firings_arrays={};
for j=1:length(varargin)
    arg0=varargin{j};
    if (isnumeric(arg0))
        if (ndims(arg0)==2)&&(size(arg0,2)>1)
            raw_arrays{end+1}=arg0;
            firings_arrays{end+1}=[];
        else
            error('Problem with argument %d',j);
        end
    elseif (iscell(arg0))
        if (length(arg0)==1)
            raw_arrays{end+1}=arg0{1};
        elseif (length(arg0)==2)
            raw_arrays{end+1}=arg0{1};
            firings_arrays{end+1}=arg0{2};
        elseif (length(arg0)==3)
            raw_arrays{end+1}=arg0{1};
            times=arg0{2}; labels=arg0{3};
            firings=zeros(3,length(times));
            firings(2,:)=times; firings(3,:)=labels;
            firings_arrays{end+1}=firings;
        end;
    elseif (isstruct(arg0))
    else
        error('Problem with argument %d',j);
    end;
end;

if (length(raw_arrays)>1)
    error('Cannot handle more than one raw arrays at this time.');
end;

if (length(raw_arrays)==1)
    mv.mode='spikespy';
    mv.raw=raw_arrays{1};
    if (length(firings_arrays{1})>1)
        mv.firings=firings_arrays{1};
    end;
    mv.sampling_freq=opts.sampling_freq;
    ms_mountainview(mv);
end;

end

function test_mv_spikespy

X=randn(4,1000);
times=randsample(size(X,2),300);
labels=ones(size(times));
opts.sampling_freq=30000;
mv_spikespy({X,times,labels},opts);

end