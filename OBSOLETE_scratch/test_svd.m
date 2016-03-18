function test_svd

X=rand(3,10);
XXt=X*X';

[U,D,V]=svd(XXt,'econ');
sqrt(D)
[U2,D2,V2]=svd(X);
D2(1:3,1:3)

U
U2

end
