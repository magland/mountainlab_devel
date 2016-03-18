function data = LoadBinary(filename,varargin)

%LoadBinary - Load data from a multiplexed binary file.
%
% jfm added dependencies right into this file - 5/22/2015
%
%  Reading a subset of the data can be done in two different manners: either
%  by specifying start time and duration (more intuitive), or by indicating
%  the position and size of the subset in terms of number of samples per
%  channel (more accurate).
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz, default = 20kHz)
%     'start'       position to start reading (in s, default = 0)
%     'duration'    duration to read (in s, default = Inf)
%     'offset'      position to start reading (in samples per channel,
%                   default = 0)
%     'samples'     number of samples (per channel) to read (default = Inf)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nChannels = 1;
precision = 'int16';
skip = 0;
frequency = 20000;
channels = [];
start = 0;
duration = Inf;
offset = 0;
nSamplesPerChannel = Inf;
time = false;
samples = false;

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+3) ' is not a property (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'start',
			start = varargin{i+1};
			if ~isdscalar(start),
				error('Incorrect value for property ''start'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if start < 0, start = 0; end
			time = true;
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>=0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			time = true;
		case 'offset',
			offset = varargin{i+1};
			if ~isiscalar(offset),
				error('Incorrect value for property ''offset'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if offset < 0, offset = 0; end
			samples = true;
		case 'samples',
			nSamplesPerChannel = varargin{i+1};
			if ~isdscalar(nSamplesPerChannel,'>=0'),
				error('Incorrect value for property ''samples'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			samples = true;
		case 'nchannels',
			nChannels = varargin{i+1};
			if ~isiscalar(nChannels,'>0'),
				error('Incorrect value for property ''nChannels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'channels',
			channels = varargin{i+1};
			if ~isivector(channels,'>=0'),
				error('Incorrect value for property ''channels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'precision',
			precision = varargin{i+1};
			if ~isstring(precision),
				error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'skip',
			skip = varargin{i+1};
			if ~isiscalar(skip,'>=0'),
				error('Incorrect value for property ''skip'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
end

% Either start+duration, or offset+size
if time && samples,
	error(['Data subset can be specified either in time or in samples, but not both (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
end

% By default, load all channels
if isempty(channels),
	channels = 1:nChannels;
end

% Check consistency between channel IDs and number of channels
if any(channels>nChannels),
	error('Cannot load specified channels (listed channel IDs inconsistent with total number of channels).');
end

% Open file
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
f = fopen(filename,'r');
if f == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end

% Size of one data point (in bytes)
sampleSize = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sampleSize = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sampleSize = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sampleSize = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sampleSize = 8;
end

% Position and number of samples (per channel) of the data subset
if time,
	dataOffset = floor(start*frequency)*nChannels*sampleSize;
	nSamplesPerChannel = floor(duration*frequency);
else
	dataOffset = offset*nChannels*sampleSize;
end

% Position file index for reading
status = fseek(f,dataOffset,'bof');
if status ~= 0,
	fclose(f);
	error('Could not start reading (possible reasons include trying to read past the end of the file).');
end

% Determine total number of samples in file
fileStart = ftell(f);
status = fseek(f,0,'eof');
if status ~= 0,
	fclose(f);
	error('Error reading the data file (possible reasons include trying to read past the end of the file).');
end
fileStop = ftell(f);
% (floor in case all channels do not have the same number of samples)
maxNSamplesPerChannel = floor((fileStop-fileStart)/nChannels/sampleSize);
frewind(f);
status = fseek(f,dataOffset,'bof');
if status ~= 0,
	fclose(f);
	error('Could not start reading (possible reasons include trying to read past the end of the file).');
end

if isinf(nSamplesPerChannel) || nSamplesPerChannel > maxNSamplesPerChannel,
	nSamplesPerChannel = maxNSamplesPerChannel;
end

% For large amounts of data, read chunk by chunk
maxSamplesPerChunk = 10000;
nSamples = nSamplesPerChannel*nChannels;
if nSamples <= maxSamplesPerChunk,
	data = LoadChunk(f,nChannels,channels,nSamples/nChannels,precision,skip);
else
	% Determine chunk duration and number of chunks
	nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
	nChunks = floor(nSamples/nSamplesPerChunk);
	% Preallocate memory
	data = zeros(nSamplesPerChannel,length(channels));
	% Read all chunks
	i = 1;
	for j = 1:nChunks,
		d = LoadChunk(f,nChannels,channels,nSamplesPerChunk/nChannels,precision,skip);
		[m,n] = size(d);
		if m == 0, break; end
		data(i:i+m-1,:) = d;
		i = i+m;
	end
	% If the data size is not a multiple of the chunk size, read the remainder
	remainder = nSamples - nChunks*nSamplesPerChunk;
	if remainder ~= 0,
		d = LoadChunk(f,nChannels,channels,remainder/nChannels,precision,skip);
		[m,n] = size(d);
		if m ~= 0,
			data(i:i+m-1,:) = d;
		end
	end
end

fclose(f);

% ---------------------------------------------------------------------------------------------------------

function data = LoadChunk(fid,nChannels,channels,nSamples,precision,skip)

if skip ~= 0,
	data = fread(fid,[nChannels nSamples],precision,skip);
else
	data = fread(fid,[nChannels nSamples],precision);
end
data=data';

if isempty(data),
	warning('No data read (trying to read past file end?)');
elseif ~isempty(channels),
	data = data(:,channels);
end

%isiscalar - Test if parameter is a scalar (integer) satisfying an optional list of tests.
%
%  USAGE
%
%    test = isiscalar(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests
%
%  EXAMPLES
%
%    % Test if x is a scalar (double)
%    isiscalar(x)
%
%    % Test if x is a strictly positive scalar (double)
%    isiscalar(x,'>0')
%
%    % Test if x is a scalar (double) included in [2,3]
%    isiscalar(x,'>=2','<=3')
%
%  NOTE
%
%    The tests ignore NaN, e.g. isiscalar(nan), isiscalar(nan,'>0') and isiscalar(nan,'<=0')
%    all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isivector, isstring,
%    islscalar, islvector, islmatrix.
%


% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isiscalar(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isiscalar">isiscalar</a>'' for details).');
end

% Test: double, scalar
test = isa(x,'double') & isscalar(x);
if ~test, return; end

% Test: integers?
test = test & round(x)==x;

% Optional tests
for i = 1:length(varargin),
	try
		if ~eval(['x' varargin{i} ';']), test = false; return; end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isiscalar">isiscalar</a>'' for details).']);
	end
end

%isivector - Test if parameter is a vector of integers satisfying an optional list of tests.
%
%  USAGE
%
%    test = isivector(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests (see examples below)
%
%  EXAMPLES
%
%    % Test if x is a vector of doubles
%    isivector(x)
%
%    % Test if x is a vector of strictly positive doubles
%    isivector(x,'>0')
%
%    % Test if x is a vector of doubles included in [2,3]
%    isivector(x,'>=2','<=3')
%
%    % Special test: test if x is a vector of doubles of length 3
%    isivector(x,'#3')
%
%    % Special test: test if x is a vector of strictly ordered doubles
%    isivector(x,'>')
%
%  NOTE
%
%    The tests ignore NaNs, e.g. isivector([500 nan]), isivector([1 nan 3],'>0') and
%    isivector([nan -7],'<=0') all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isiscalar, isstring,
%    islscalar, islvector, islmatrix.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isivector(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isivector">isivector</a>'' for details).');
end

% Test: double, vector
test = isa(x,'double') & isvector(x);

% Ignore NaNs
x = x(~isnan(x));

% Test: integers?
test = test & all(round(x)==x);

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if length(x) ~= str2num(varargin{i}(2:end)), test = false; return; end
		elseif isstring(varargin{i},'>','>=','<','<='),
			dx = diff(x);
			if ~eval(['all(0' varargin{i} 'dx);']), test = false; return; end
		else
			if ~eval(['all(x' varargin{i} ');']), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isivector">isivector</a>'' for details).']);
	end
end

%isstring - Test if parameter is an (admissible) character string.
%
%  USAGE
%
%    test = isstring(x,string1,string2,...)
%
%    x              item to test
%    string1...     optional list of admissible strings
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isivector, isiscalar.
%

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isstring(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isstring">isstring</a>'' for details).');
end

test = true;

if ~ischar(x),
	test = false;
	return;
end

if isempty(varargin), return; end

for i = 1:length(varargin),
	if strcmp(x,varargin{i}), return; end
end

test = false;
