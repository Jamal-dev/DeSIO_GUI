% =================================================================================================================
function DeSiO_Aero_plot_qw(currDir,strfilename,wake,time)
cd = currDir;

model = uvlm_readmodel();
t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);

for i = 1:length(time)
    [m,iz]     = min(abs(t-time(i)));
    inzstep(i) = iz;
    inz(i)     = (t(2)-t(1))*t(1) - 1 + iz;
end

aa = 1; coord =[];
for k = 1:size(wake,2)
    cell_qw = read_wake([model.strSimName '_uvlm_qw' num2str(wake(k).nr) '.dres']);
    npwe    = wake(k).node_per_wake_edge;
    nodes   = wake(k).nodes;
    for i = 1:length(nodes)
        for j = 1:length(time)
            node_i = npwe*(inz(j))+nodes(i);
            inzPos = 3*(node_i-1)+1:3*(node_i-1)+3;
            coord(j,3*(aa-1)+1:3*(aa-1)+3) = cell_qw{inzstep(j)}(inzPos);
        end
        aa = aa + 1;
    end
end
fid = fopen([strfilename '_wake.txt'],'w');
for i = 1:size(coord,1)
    for j = 1:size(coord,2)
        if j~=size(coord,2)
            fprintf(fid,'%10.8f\t',coord(i,j));
        else
            fprintf(fid,'%10.8f\n',coord(i,j));
        end
    end
end
fclose(fid);

function qw = read_wake(strfilename)
    qw = {};
    fid = fopen(strfilename,'r');
    %tline = fgetl(fid);
    while ~feof(fid)
        tline = fgetl(fid);
        qw{end+1} = sscanf(tline, '%f');
    end
    fclose(fid);
return



