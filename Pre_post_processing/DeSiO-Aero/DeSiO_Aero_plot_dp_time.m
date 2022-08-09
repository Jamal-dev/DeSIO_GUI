% =================================================================================================================
function DeSiO_Aero_plot_dp_time(currDir,strfilename,rings,time,rgbColor,plotStyle)
cd = currDir;

% check, if input exists
if ~exist('rgbColor','var')
    rgbColor = abs(rand(length(rings),3));
end

if ~exist('plotStyle','var')
    for i = 1:length(rings)
        plotStyle{i} = '-';
    end
end

% plot lift coeficient
model = uvlm_readmodel();
t     = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
dp    = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);

field = dp;

for i = 1:length(time)
    [m,iz] = min(abs(t-time(i)));
    inz(i)   = iz;
end

figure(); hold on; grid on; hold on;
set(gca,'fontweight','bold','fontsize',12);
xlabel('time s'); ylabel(['\Deltap']);

leg = {};
for i = 1:length(rings)
    plot(t(:),field(:,rings(i)),'linewidth',1.5,'color',rgbColor(i,:),'linestyle',plotStyle{i});
    leg(end+1) = {['BE ' num2str(rings(i))]};
end
legend(leg,'Location','southoutside','Orientation','horizontal');

for j = 1:length(rings)
    for i = 1:length(inz)
        plot(t(inz(i)),field(inz(i),rings(j)),'.','color',rgbColor(j,:),'markersize',14);
    end
end
print([strfilename '_dp_vs_time'],'-dpng', '-r500');
close all;

fid = fopen([strfilename '_dp.txt'],'w');
for i = 1:length(inz)
    for j = 1:length(rings)
        if j~=length(rings)
            fprintf(fid,'%10.8f\t',field(inz(i),rings(j)));
        else
            fprintf(fid,'%10.8f\n',field(inz(i),rings(j)));
        end
    end
end
fclose(fid);
return