function test_get_pca_features

M=4; N=10;
Y=zeros(M,N);
Y(1,:)=1;
Y(2,:)=(0:N-1);
Y(3,:)=(0:N-1).^2;
Y(4,:)=(0:N-1).^3;
Y*Y'
Y=reshape(Y,1,M,N);
ms_event_features(Y,4)

end