% Test consistency between ms_whiten and mscmd_whiten

X=rand(4,100);
Y1=ms_whiten(X);

writemda(X,'tmp_X.mda');
mscmd_whiten('tmp_X.mda','tmp_X_white.mda');
Y2=readmda('tmp_X_white.mda');

disp(max(abs(Y1(:)-Y2(:))));