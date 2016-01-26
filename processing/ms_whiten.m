function X=ms_whiten(X)

mu = mean(X,2); 
X = X-repmat(mu,1,size(X,2));
[U,D,V] = svd(X*X');
X=sqrt(size(X,2)-1)*U*sqrt(inv(D))*(U'*X);

end
