function compare_rotor_velocity()
% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
close all
format short

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

currDir       = cd;
path_desio    = ['\DeSiO\15MW_pitch10_vel12_nspan3\DeSiO'];
path_ofast_BD = ['\OpenFAST_BeamDyn'];
path_ofast_ED = ['\OpenFAST_ElastoDyn'];

str_ofast_BD_filename = 'IEA-15-240-RWT-Onshore';
str_ofast_ED_filename = 'IEA-15-240-RWT-Onshore';

ndeltafactor = 1/10;

lineSize   = 1;
markerSize = 6;
FontSize   = 12;

flag_plot_time_hist = 1;
flag_plot_mean_val  = 1;

flag_plot_De        = 1;
flag_plot_ED        = 1;
flag_plot_BD        = 1;

flag_save_fig       = 0;

flag_delete_De      = 0;
flag_delete_ED      = 1;
flag_delete_BD      = 0;

var_name    = 'rotorspeed';
var_name_ED = 'RtSpeed';
var_name_BD = 'RtSpeed';

nd_ax = -[-7.03233176e-01;7.03233176e-01;1.04528463e-01]; % axis in local director cos of node
node  = 2;

% column for results of wing tip displacement
arrcolumn = [];
nplotDe  = 0; nplotBD = 0; nplotED = 0; 

%% read DeSiO-Results... tower top
if flag_plot_De == 1
    cd([currDir '' path_desio]);
    % read model
    model_fsi       = fsi_readmodel;
    model_structure = structure_readmodel;
    if isempty(model_fsi) && not(isempty(model_structure))
        model = model_structure;
    elseif isempty(model_structure) && not(isempty(model_fsi))
        model = model_fsi;
    else
        disp(['No model found!']);
        return;
    end
    time = load([model.strSimName '_t.dres']); time = time(:,1);
    if exist(['res_' var_name '_de.mat'], 'file') && flag_delete_De == 0
        res_de = load(['res_' var_name '_de'],'-mat');
        res_de = res_de.res_de;
    else
        q       = load([model.strSimName '_q.dres']);
        s       = load([model.strSimName '_v.dres']);
        res_de  = fun_extract_rotor_speed(q,s,time,node,nd_ax);
        save(['res_' var_name '_de'],'res_de');
    end
    nplotDe = 1;
    inztime2plotDe = [1:1/ndeltafactor:size(time,1)];
end

%% read OpenFAST ElastoDyn results
if flag_plot_ED == 1
    cd([currDir '' path_ofast_ED]);
    if exist(['res_' var_name  '_ofastED.mat'], 'file') && flag_delete_ED == 0
        ofastED_res = load(['res_' var_name  '_ofastED'],'-mat'); ofastED_res = ofastED_res.ofastED_res;
        strInpVarED = load('strInpVarED','-mat'); 
        strInpVarED = strInpVarED.strInpVarED;
    else
        [ofastED_res,strInpVarED] = fun_LoadOpenFastResFile([str_ofast_ED_filename '.out']);
        save(['res_' var_name  '_ofastED.mat'],'ofastED_res');
        save('strInpVarED.mat','strInpVarED');
    end

    for i = 1:size(strInpVarED{1},1)
        if strncmp(var_name_ED,strInpVarED{1}(i),9) == 1; arrcolumnED(1,1) = i; end
    end
    cd(currDir);
    nplotED = size(arrcolumnED,1);
    inztime2plotED = [1:1/ndeltafactor:size(ofastED_res,1)];
end

%% read OpenFAST BeamDyn results
if flag_plot_BD == 1
    cd([currDir '' path_ofast_BD]);
    if exist(['res_' var_name  '_ofastBD.mat'], 'file') && flag_delete_BD == 0
        ofastBD_res = load(['res_' var_name  '_ofastBD'],'-mat'); ofastBD_res = ofastBD_res.ofastBD_res;
        strInpVarBD = load('strInpVarBD','-mat'); 
        strInpVarBD = strInpVarBD.strInpVarBD;
    else
        [ofastBD_res,strInpVarBD] = fun_LoadOpenFastResFile([str_ofast_BD_filename '.out']);
        save(['res_' var_name  '_ofastBD.mat'],'ofastBD_res','-v7.3');
        save('strInpVarBD.mat','strInpVarBD');
    end

    for i = 1:size(strInpVarBD{1},1)
        if strncmp(var_name_BD,strInpVarBD{1}(i),size(var_name_BD,2)) == 1; arrcolumnBD(1,1) = i; end
    end
    cd(currDir);
    nplotBD = size(arrcolumnBD,1);
    inztime2plotBD = [1:1/ndeltafactor:size(ofastBD_res,1)];
end

if flag_plot_BD == 1 && flag_plot_ED == 0 && flag_plot_De == 0; nsize = min([nplotBD]); end
if flag_plot_BD == 0 && flag_plot_ED == 1 && flag_plot_De == 0; nsize = min([nplotED]); end
if flag_plot_BD == 0 && flag_plot_ED == 0 && flag_plot_De == 1; nsize = min([nplotDe]); end
if flag_plot_BD == 1 && flag_plot_ED == 1 && flag_plot_De == 0; nsize = min([nplotBD,nplotED]); end
if flag_plot_BD == 1 && flag_plot_De == 1 && flag_plot_ED == 0; nsize = min([nplotBD,nplotDe]); end
if flag_plot_ED == 1 && flag_plot_De == 1 && flag_plot_BD == 0; nsize = min([nplotED,nplotDe]); end
if flag_plot_ED == 1 && flag_plot_De == 1 && flag_plot_BD == 1; nsize = min([nplotBD,nplotED,nplotDe]); end

% Function to extract dof-solution for node from solution.m file
for i = 1:nsize

    if flag_plot_time_hist == 1
        
%% Plotting time series
        strWhatToComp = []; strleg = {};
        handlesPlot1 = []; handlesPlot2 = []; handlesPlot3 = [];

        figure(); hold on; grid on;
        hold on; grid on; title([var_name]); ylabel('-'); xlabel('time in s');

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];
            plot(time,res_de(i).val(:,1),'-r','linewidth',lineSize);
            handlePlotDe1 = plot(time(inztime2plotDe),res_de(i).val(inztime2plotDe,1),'.','color','k','marker','o','markersize',markerSize);
            handlesPlot1 = [handlesPlot1,handlePlotDe1]; 
            strleg = [strleg,{'DeSiO'}];
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];
            plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,1)),'-r','linewidth',lineSize);
            handlePlotED1 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,1)),'.','color','k','marker','+','markersize',markerSize);
            handlesPlot1 = [handlesPlot1,handlePlotED1]; 
            strleg = [strleg,{'OpenFAST ED'}];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];
            plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,1)),'-r','linewidth',lineSize);
            handlePlotBD1 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,1)),'.','color','k','marker','s','markersize',markerSize);
            handlesPlot1 = [handlesPlot1,handlePlotBD1]; 
            strleg = [strleg,{'OpenFAST BD'}];
        end

        if isempty(strleg) ~= 1
            leg1 = legend(handlesPlot1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
        end

        if flag_save_fig == 1
            fig1 = gcf; 
            frame_h = get(handle(fig1),'JavaFrame'); set(frame_h,'Maximized',1);
            saveFigure(gcf,['compare_' var_name '_' strWhatToComp '.fig']);
            print(['compare_' var_name '_' strWhatToComp],'-dpng', '-r500');
            set(frame_h,'Maximized',0);
        end
    end
    
    if flag_plot_mean_val == 1
%% plotting mean value
        strleg = {}; strWhatToComp = []; arrColor = [];
        val_m = [];

        figure(); hold on; grid on; title([var_name]); ylabel('mean value'); 

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];
            val_m_de = mean(res_de(i).val(:,1)); 
            val_m    = [val_m,val_m_de];
            arrColor = [arrColor;[1,0,0]]; strleg = [strleg,{'DeSiO'}];
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];
            val_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,1)));
            val_m = [val_m,val_m_ofastED];
            arrColor = [arrColor;[0,1,0]]; strleg = [strleg,{'OpenFAST ED'}];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];
            val_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,1)));
            val_m = [val_m,val_m_ofastBD];
            arrColor = [arrColor;[0,0,1]]; strleg = [strleg,{'OpenFAST BD'}];
        end

        % bar plot for mean values
        if isempty(val_m) ~= 1; bar1 = bar(vertcat([1:length(val_m)],nan(length(val_m))),vertcat(val_m,nan(length(val_m)))); end
        if isempty(strleg) ~= 1
            leg1 = legend(bar1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
        end
        
        if flag_save_fig == 1
            fig2 = gcf; 
            frame_h = get(handle(fig2),'JavaFrame'); set(frame_h,'Maximized',1);
            saveFigure(gcf,['compare_' var_name '_mean_' strWhatToComp '.fig']);
            print(['compare_' var_name '_mean_' strWhatToComp],'-dpng', '-r500');
            set(frame_h,'Maximized',0);
        end
    end
end

return