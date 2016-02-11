function test_t_mixture

nu=5;
n=100;

mu1=zeros(1,n);
mu2=zeros(1,n); mu2(1)=0.5;
sigma=1;

p_s=[];

for ct=1:50
    v=randn(1,n)*12;
    dist1=sqrt(sum((v-mu1).^2));
    dist2=sqrt(sum((v-mu2).^2));
    %   exp(-0.5*dist1^2/sigma^2)/(exp(-0.5*dist1^2/sigma^2)+exp(-0.5*dist2^2/sigma^2))
    % = exp( -0.5*dist1^2/sigma^2 - log(exp(-0.5*dist1^2/sigma^2)+exp(-0.5*dist2^2/sigma^2)) )

    lognumerator=-0.5*dist1.^2/sigma^2
    logdenominator=logsumexp(-0.5*[dist1,dist2].^2/sigma^2);

    p=exp(lognumerator-logdenominator)
    p_s(end+1)=p;
end;

figure; hist(p_s,0:0.05:1);

end

