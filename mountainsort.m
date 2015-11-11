function mountainsort(opts)

%If there are no inputs, we run the test
if (nargin==0) mountainsort_example1; return; end;

total_timer=tic;

%Setup the environment, create the working directories, etc
mountainsort_setup(opts);

%Step 1: Preprocess
fprintf('\n*** Step 1: Preprocess...\n');
data=step1_preprocess(opts);
fprintf('Writing raw.mda...\n');
writemda(data.X,[opts.output_path,'/raw.mda']);

%Step 1a: Prepare
fprintf('\n*** Step 1a: Prepare...\n');
step1a_prepare(opts);

%Step 2: Detect
fprintf('\n*** Step 2: Detect...\n');
step2_detect(opts,data);

%Step 3: Cluster
fprintf('\n*** Step 3: Cluster...\n');
step3_cluster(opts,data);

%Step 4: Consolidate
fprintf('\n*** Step 4: Consolidate...\n');
[times,labels]=step4_consolidate(opts);

%Cross-correlogram
fprintf('\n*** Computing cross-correlograms ***\n');
max_dt=1500;
CC=compute_cross_correlograms(times,labels,max_dt);
path0=[opts.output_path,'/cross-correlograms.mda'];
fprintf('Writing %s...\n',path0);
writemda(cross_correlograms_to_mda(CC),path0);

fprintf('Total processing time: %.2f sec\n',toc(total_timer));

end

function mountainsort_setup(opts)

basepath=fileparts(mfilename('fullpath'));

fprintf('Setting up... basepath=%s...\n',basepath);
addpath([basepath,'/view']);
addpath([basepath,'/util']);
addpath([basepath,'/processing']);
addpath([basepath,'/isosplit']);

fprintf('Verifying input parameters in opts...\n');
if (~isfield(opts,'timepoints')) opts.timepoints=[]; end;
required_fields={...
    'working_path',...
    'channels',...
    'raw_dat',...
    'raw_mda',...
};
for j=1:length(required_fields)
    if (~isfield(opts,required_fields{j}))
        error('opts.%s is a required input parameter.',required_fields{j});
    end;
end;

fprintf('Creating working directories in %s...\n',opts.working_path);
if (~exist(opts.working_path,'dir')) mkdir(opts.working_path); end;
if (~exist([opts.working_path,'/detect'],'dir')) mkdir([opts.working_path,'/detect']); end;
if (~exist([opts.working_path,'/cluster'],'dir')) mkdir([opts.working_path,'/cluster']); end;
if (~exist(opts.output_path,'dir')) mkdir(opts.output_path); end;

end

function ret=cross_correlograms_to_mda(CC)
K=size(CC,1);

ct=0;
for k1=1:K
for k2=1:K
ct=ct+length(CC{k1,k2});
end;
end;
ret=zeros(3,ct);

ct=0;
for k1=1:K
for k2=1:K
ret(1,ct+1:ct+length(CC{k1,k2}))=k1;
ret(2,ct+1:ct+length(CC{k1,k2}))=k2;
ret(3,ct+1:ct+length(CC{k1,k2}))=CC{k1,k2};
ct=ct+length(CC{k1,k2});
end;
end;

end

