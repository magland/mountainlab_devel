function print_hdf5_info(path)

fprintf('%s\n',path);
info=hdf5info(path);
print_hdf5_group(info.GroupHierarchy,path);

end


function print_hdf5_group(G,path)
for j=1:length(G.Groups)
    print_hdf5_group(G.Groups(j),path);
end;
for j=1:length(G.Datasets)
    print_hdf5_dataset(G.Datasets(j),path);
end;
end

function print_hdf5_dataset(D,path)
dtype=D.Datatype(1).Class;
fprintf('%s: ',D.Name);
if (strcmp(dtype,'H5T_STRING'));
    tmp=hdf5read(path,D.Name);
    str=tmp.Data;
    fprintf('"%s"',str);
    fprintf('\n');
else
    if (length(D.Dims)==0)
        tmp=hdf5read(path,D.Name);
        fprintf('%g\n',tmp);
    else
        fprintf('%s %s\n',dtype,num2str(D.Dims));
    end;
end;
end
