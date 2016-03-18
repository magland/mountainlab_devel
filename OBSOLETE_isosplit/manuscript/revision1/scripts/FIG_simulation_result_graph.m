function FIG_simulation_result_graph

close all;

mfile_path=fileparts(mfilename('fullpath'));

str=fileread(sprintf('%s/../isosplit-revision1.tex',mfile_path));
lines=strsplit(str,'\n');
DATA1=zeros(0,3);
DATA2=zeros(0,3);
for j=1:length(lines)
	line=lines{j};
	if (length(strfind(line,'\pm'))>0)&&(length(strfind(line,' & '))>2)
		disp(line);
		line=strrep(line,'$','');
		line=strrep(line,' ','');
		line=strrep(line,'\%','');
		line=strrep(line,'\\','');
		strings=strsplit(line,'&');
		tmp1=[];
		tmp2=[];
		for k=2:length(strings)
			str=strings{k};
			vals=strsplit(str,'\\pm');
			tmp1(end+1)=eval(vals{1});
			tmp2(end+1)=eval(vals{2});
		end;
		DATA1(end+1,:)=tmp1;
		DATA2(end+1,:)=tmp2;
	end;
end;

DATA1=reshape(DATA1,5,size(DATA1,1)/5,size(DATA1,2));
DATA2=reshape(DATA2,5,size(DATA2,1)/5,size(DATA2,2));

for cc=1:3
	figure('name',sprintf('cc=%d',cc)); set(gcf,'Color','w');
	barweb(DATA1(:,:,cc)',DATA2(:,:,cc)'); hold on;
	ylabel('Avg. Accuracy (%)');
	legend('ISO-SPLIT','K-means (K known)','GMM (K known)','GMM (BIC)','DBSCAN','location','northeastoutside');
	set(gca,'xticklabel',{'S1: Isotropic','S2: Anisotropic','S3: Non-Gaussian','S4: Packed','S5: High-dimensional'});
	set(gcf,'PaperPosition',[0,0,12,5]);
	if (cc==2)
		print(sprintf('%s/../simulation_result_graph.eps',mfile_path),'-depsc2');
	end;
end;

end
