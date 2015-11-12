function [TIMES,LABELS]=step5_consolidate(opts)

timerA=tic;

cluster_times_prefix=[opts.working_path,'/cluster/times_'];
cluster_labels_prefix=[opts.working_path,'/cluster/labels_'];
cluster_templates_prefix=[opts.working_path,'/cluster/templates_'];
consolidate_times_path=[opts.output_path,'/times.mda'];
consolidate_labels_path=[opts.output_path,'/labels.mda'];
consolidate_templates_path=[opts.output_path,'/templates.mda'];
consolidate_primary_channels_path=[opts.output_path,'/primary_channels.mda'];
adjacency_path=[opts.output_path,'/adjacency.mda'];

AM=readmda(adjacency_path);

M=size(AM,1);

TIMES=[];
LABELS=[];
TEMPLATES=zeros(M,opts.clip_size,0);
PRIMARY_CHANNELS=[];

current_label=1;

for j=1:M
    fname_cluster2_times=[cluster_times_prefix,sprintf('%d.mda',j)];
    fname_cluster2_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster2_templates=[cluster_templates_prefix,sprintf('%d.mda',j)];
    times=readmda(fname_cluster2_times);
    labels=readmda(fname_cluster2_labels);
    WF=readmda(fname_cluster2_templates);
    
    if (size(WF(:)>1))
        sizes=squeeze(sum(WF.^2,2));
        max_sizes=(max(sizes,[],1));
        rel_sizes=sizes(j,:)./max_sizes;
        labels_to_use=find(rel_sizes>=0.9);

        for ii=1:length(labels_to_use)
            k=labels_to_use(ii);
            time_indices=find(labels==k);
            TIMES=[TIMES,times(time_indices)];
            LABELS=[LABELS,ones(size(labels(time_indices)))*current_label];
            TEMPLATES=cat(3,TEMPLATES,WF(:,:,k));
            PRIMARY_CHANNELS=[PRIMARY_CHANNELS,j];
            current_label=current_label+1;
        end;
    end;
end;

[TIMES,sort_inds]=sort(TIMES);
LABELS=LABELS(sort_inds);

fprintf('Writing %s...\n',consolidate_times_path);
writemda(TIMES,consolidate_times_path);
fprintf('Writing %s...\n',consolidate_labels_path);
writemda(LABELS,consolidate_labels_path);
fprintf('Writing %s...\n',consolidate_templates_path);
writemda(TEMPLATES,consolidate_templates_path);
fprintf('Writing %s...\n',consolidate_primary_channels_path);
writemda(PRIMARY_CHANNELS,consolidate_primary_channels_path);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end
