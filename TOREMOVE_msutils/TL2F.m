function firings=TL2F(times,labels)
firings=zeros(3,length(times));
firings(2,:)=times;
firings(3,:)=labels;
end