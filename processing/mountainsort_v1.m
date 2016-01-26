function mountainsort_v1(rawfile_path,outputdir_path,opts)

timerA=tic;

path0=outputdir_path;

if nargin<3, opts=struct; end;

if (~isfield(opts,'adjacency'))
    error('opts.adjacency not found.');
end;
writemda_if_not_equal(opts.adjacency,[path0,'/adjacency.mda']);

if isfield(opts,'o_filter')
    mscmd_bandpass_filter(rawfile_path,[path0,'/pre1.mda'],opts.o_filter);
end;
if isfield(opts,'o_whiten')
    %mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],opts.o_whiten);
    writemda(ms_whiten(readmda([path0,'/pre1.mda'])),[path0,'/pre2.mda']);
end;
if isfield(opts,'o_filter2')
    mscmd_bandpass_filter([path0,'/pre2.mda'],[path0,'/pre3.mda'],opts.o_filter2);
end;
if isfield(opts,'o_detect')
    if ~isfield(opts.o_detect,'use_pre1')
        opts.o_detect.use_pre1=0;
    end;
end;
if isfield(opts,'o_split_clusters')
    if ~isfield(opts.o_split_clusters,'use_pre1')
        opts.o_split_clusters.use_pre1=0;
    end;
end;
if ~isfield(opts.o_detect,'use_pre1')
    opts.o_split_clusters.use_pre1=0;
end;

if opts.o_detect.use_pre1
    mscmd_detect([path0,'/pre1.mda'],[path0,'/detect.mda'],opts.o_detect);
else
    mscmd_detect([path0,'/pre3.mda'],[path0,'/detect.mda'],opts.o_detect);
end;
mscmd_features([path0,'/pre3.mda'],[path0,'/detect.mda'],[path0,'/adjacency.mda'],[path0,'/features.mda'],opts.o_features);
mscmd_cluster([path0,'/features.mda'],[path0,'/clusters1.mda'],opts.o_cluster);

if (isfield(opts,'o_split_clusters'))
    if opts.o_split_clusters.use_pre1
        mscmd_split_clusters([path0,'/pre1.mda'],[path0,'/clusters1.mda'],[path0,'/clusters2.mda'],opts.o_split_clusters);
    else
        mscmd_split_clusters([path0,'/pre3.mda'],[path0,'/clusters1.mda'],[path0,'/clusters2.mda'],opts.o_split_clusters);
    end;
else
    mscmd_copy([path0,'/clusters1.mda'],[path0,'/clusters2.mda']);
end;
mscmd_templates([path0,'/pre3.mda'],[path0,'/clusters2.mda'],[path0,'/templates1.mda'],opts.o_templates);
if isfield(opts,'o_consolidate')
mscmd_consolidate([path0,'/clusters2.mda'],[path0,'/templates1.mda'],[path0,'/clusters3.mda'],[path0,'/templates2.mda'],[path0,'/load_channels.mda'],opts.o_consolidate);
mscmd_fit([path0,'/pre3.mda'],[path0,'/clusters3.mda'],[path0,'/templates2.mda'],[path0,'/clusters.mda'],opts.o_fit);
else
mscmd_fit([path0,'/pre3.mda'],[path0,'/clusters2.mda'],[path0,'/templates1.mda'],[path0,'/clusters.mda'],opts.o_fit);    
end;



if (~isempty(rawfile_path))
    mscmd_templates(rawfile_path,[path0,'/clusters.mda'],[path0,'/templates_raw.mda'],opts.o_templates);
end;
mscmd_templates([path0,'/pre3.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],opts.o_templates);
mscmd_templates([path0,'/pre1.mda'],[path0,'/clusters.mda'],[path0,'/templates_pre1.mda'],opts.o_templates);

mscmd_extract_clips([path0,'/pre3.mda'],[path0,'/clusters.mda'],[path0,'/clips.mda'],[path0,'/clips_index.mda'],opts.o_extract_clips);
mscmd_extract_clips([path0,'/pre1.mda'],[path0,'/clusters.mda'],[path0,'/clips_pre1.mda'],[path0,'/clips_index.mda'],opts.o_extract_clips);

mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross_correlograms.mda'],opts.o_cross_correlograms.max_dt);

fprintf('Elapsed time: %g sec\n',toc(timerA));

end

function writemda_if_not_equal(X,fname)
if (~mda_equals_file(X,fname)) writemda(X,fname); end;
end

function ret=mda_equals_file(X,fname)
if (exist(fname,'file'))
    Y=readmda(fname);
    ret=isequal(X,Y);
else
    ret=0;
end;
end

