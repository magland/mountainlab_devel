function mv_view_clusters(X,labels)
%MV_VIEW_CLUSTERS - 3D view of clusters in feature space, colored by
%labels
%
% Syntax:  mv_view_clusters(features,labels)
%
% Inputs:
%    features - RxL array of features (R>=3), e.g., optained by
%               mscmd_features or ms_event_features
%    labels - 1xL array of integer labels (controls the colors of the data
%    points)
%
% Other m-files required: ms_mountainsort
%
% See also: ms_view_clusters, ms_event_features

% Author: Jeremy Magland
% Mar 2016. Last revision: 9 Mar 2016

if nargin<1, test_mv_view_clusters; return; end;

mv.mode='view_clusters';
mv.data=X;
if (nargin>=2)
    mv.labels=labels;
end;
ms_mountainview(mv);

end

function test_mv_view_clusters

X=zeros(3,0);
labels=zeros(1,0);
X=cat(2,X,randn(3,5000)); labels=cat(2,labels,ones(1,5000));
X=cat(2,X,randn(3,5000)+3); labels=cat(2,labels,ones(1,5000)*2);

mv_view_clusters(X,labels);
ms_view_clusters(X,labels);

end