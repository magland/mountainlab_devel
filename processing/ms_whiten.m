function X=ms_whiten(X)

if nargin<1, test_ms_whiten; return; end;

[M,N]=size(X);

mu = mean(X,2); 
X = X-repmat(mu,1,size(X,2));
[U,D,V] = svd(X);
D(D~=0)=1./D(D~=0);
X=U*D(1:M,1:M)*(U'*X);

end

function X=ms_whiten_XXt(X)
%used to show that we get the same result

mu = mean(X,2); 
X = X-repmat(mu,1,size(X,2));
[U,D,V] = svd(X*X');
X=U*sqrt(inv(D))*(U'*X);

end

function X=ms_whiten_mscmd(X)
writemda(X,'X_tmp.mda');
mscmd_whiten('X_tmp.mda','X_tmp_white.mda');
X=readmda('X_tmp_white.mda');
end

function test_ms_whiten

X=rand(4,1000);
Y1=ms_whiten(X);
Y1*Y1'
Y2=ms_whiten_XXt(X);
Y2*Y2'
Y3=ms_whiten_mscmd(X);
Y3*Y3'

max(abs(Y1(:)-Y2(:)))
max(abs(Y1(:)-Y3(:)))

end