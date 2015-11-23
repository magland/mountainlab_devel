function test_fast_pca

close all;

M=1800; N=10000; npca=6;
X=randn(M,N)*0.1;
%X=X+ones(size(X))*0.1;
X(:,1:2:end)=X(:,1:2:end)+repmat((1:M)'/M,1,length(1:2:N))*0.1;
X(:,2:2:end)=X(:,2:2:end)+repmat((M:-1:1)'/M,1,length(2:2:N))*0.1;

X0=X;
lambdas=[];
for k=1:npca
disp('.');
[v,lambda]=get_top_comp(X0);
X0=subtract_component(X0,v);
lambdas(end+1)=lambda;
end

tic;
[U,D]=eig(X*X');
[d,I]=sort(diag(D),'descend');
U=U(:,I);
z = U'*X;
toc

plot(1:npca,lambdas,'b',1:npca,d(1:npca),'r',1:npca,lambdas'-d(1:npca),'k');

end

function X=subtract_component(X,v)
[M,N]=size(X);
for n=1:N
    X(:,n)=X(:,n)-(X(:,n)'*v)*v;
end;
end

function [v,lambda]=get_top_comp(X)

%Y=X*X';

M=size(X,1);
v=randn(M,1);
maxit=100;
numit=0;
while 1
    v_last=v;
    v=X*(X'*v);
    %v=Y*v;
    norm_v=sqrt(v'*v);
    v=v/norm_v;
    diff=sqrt(sum((v-v_last).^2));
    %if (diff<0.001) break; end;
    numit=numit+1;
    if (numit>=maxit) break; end;
end;
lambda=norm_v;
end


%function [z U] = features_pca(X)
% FEATURES_PCA - get principal component analysis feature vectors, for Ns large
% Jeremy's version using X X^T, faster for large Ns, but limits to 8 digits?
% Assumes that M*Nt < Ns otherwise it's slower than plain SVD on X.
% Tries to standardize signs of z.
% todo: investigate QR or LQ for highly fat matrix SVD case.

% [M Nt Ns] = size(X);           % Get some dimensions
% MM=M*Nt; X=reshape(X,MM,Ns);   % collapse channel and time dimensions
% [U,D] = eig(X*X');   % takes O(MM^2*(MM+Ns)). Note eig faster than svd.
% [d,I] = sort(diag(D),'descend'); U = U(:,I);  % sort eigenvectors
% U = bsxfun(@times, U, std_signs(U));   % std signs of col vecs
% z = U'*X;   % get all components in O(MM^2*Ns).   sing vals = sqrt(diag(D))
% % sqrt(d(1:10)), U(1:10,1)   % few singular values & 1st left vec