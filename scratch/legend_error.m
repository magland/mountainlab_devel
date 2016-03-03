function legend_error

figure;
for j=1:5
    x=rand(1,10);
    y=rand(1,10);
    z=rand(1,10);
    scatter3(x,y,z); hold on;
end;

%MATLAB bug: Must use two output parameters for this to work - thanks Alex
[~,~]=legend('A','B','C','D','E');

