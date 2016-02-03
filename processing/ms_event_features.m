function [z subspace] = ms_event_features(X,num_features,o)

if nargin<3, o = []; end
if ~isfield(o,'fmethod'), o.fmethod='pca'; end  % default
[M T Ns] = size(X);
tic;
% call requested alg...
if strcmp(o.fmethod,'pca'), [z U] = features_pca(X);
  subspace = reshape(U,[size(X,1),size(X,2),size(U,2)]);
elseif strcmp(o.fmethod,'pcachan')
  z = zeros(M*T,Ns); U = zeros(M*T,M*T);
  for m=1:M                % PCA for each channel separately
    r = m+(0:T-1)*M;  % row indices to write to: all 1st dims together, etc.
    [z(r,:) U(r,r)]= features_pca(X(m,:,:));
  end
  subspace = reshape(U,[size(X,1),size(X,2),size(U,2)]);
elseif strcmp(o.fmethod,'cen'), z = features_spatial_centroid(X,o.d);
elseif strcmp(o.fmethod,'raw')
  z = reshape(permute(X,[2 1 3]),[M*T Ns]); % matrix: events are cols
else, error('unknown fmethod in features!');
end
%fprintf('features done in %.3g s\n',toc)
%%%%%

z=z(1:num_features,:);
subspace=subspace(:,:,1:num_features);

function [z U] = features_pca(X)
% FEATURES_PCA - get principal component analysis feature vectors, for Ns large
% Jeremy's version using X X^T, faster for large Ns, but limits to 8 digits?
% Assumes that M*Nt < Ns otherwise it's slower than plain SVD on X.
% Tries to standardize signs of z.
% todo: investigate QR or LQ for highly fat matrix SVD case.

[M Nt Ns] = size(X);           % Get some dimensions
MM=M*Nt; X=reshape(X,MM,Ns);   % collapse channel and time dimensions
[U,D] = eig(X*X');   % takes O(MM^2*(MM+Ns)). Note eig faster than svd.
[d,I] = sort(diag(D),'descend'); U = U(:,I);  % sort eigenvectors
U = bsxfun(@times, U, std_signs(U));   % std signs of col vecs
z = U'*X;   % get all components in O(MM^2*Ns).   sing vals = sqrt(diag(D))
% sqrt(d(1:10)), U(1:10,1)   % few singular values & 1st left vec

function z = features_pca_crude(X)
% FEATURES_PCA_CRUDE - get principal component analysis feature vectors.
% Alex's version, plain SVD on A, slower for large Ns, full e_mach.
% Tries to standardize signs of z.
[M T Ns] = size(X);
X = reshape(permute(X,[2 1 3]),[M*T Ns]); % matrix: events are cols
% (the permute has no effect for PCA but is my standard way to unfold to matrix)
[U S V] = svd(X,'econ');     % O((M*T)^2 * Ns) ? not sure for a fat matrix
S = bsxfun(@times, S, std_signs(U));   % std signs of col vecs
z = S*V';  % silly, could use repmat to make O(M*T*Ns)

function s = std_signs(U)
% standardized signs from col vecs of U. Barnett 12/23/14
[m n] = size(U);
s = nan(1,n);
for j=1:n
  [~,i] = max(abs(U(:,j)));   % index of max entry in jth col vec of U
  s(j) = sign(U(i,j));
end

function z = features_spatial_centroid(X,d)
% electrode location features, idea from Prentice et al 2011. Barnett 12/22/14
% Inputs:
%  X - usual 3D single-event array
%  d - raw EC data struct, for electrode loc info
% Outputs:
%  z - feature vectors as 3*Ns array
[M T Ns] = size(X);
z = nan(3,Ns);  % allocate output
for i=1:Ns         % loop over events...
  w = max(-X(:,:,i),[],2);  % weight by negative voltage peak - Prentice
  z(3,i) = max(w);          % 3rd coord is max peak height
  w = w(:)/sum(w);        % normalized weights col vec
  z(1:2,i) = d.electrodelocs * w;  % weighted spatial electrode locs in xy
end


