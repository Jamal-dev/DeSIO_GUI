function compare_wt_tip_deformation()
% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
close all
format short

addpath(genpath('U:\PhD-2021-12_Christian_Hente\Pre_post_processing'));

currDir       = cd;
path_desio    = ['\DeSiO\15MW_pitch10_vel12_nspan3\DeSiO'];
path_ofast_BD = ['\OpenFAST_BeamDyn'];
path_ofast_ED = ['\OpenFAST_ElastoDyn'];

str_ofast_BD_filename = 'IEA-15-240-RWT-Onshore';
str_ofast_ED_filename = 'IEA-15-240-RWT-Onshore';

var_name    = 'tip_displacement';

ndeltafactor = 1/5;

lineSize   = 1;
markerSize = 6;
FontSize   = 12;

flag_plot_time_hist = 1;
flag_plot_mean_val  = 1;

flag_plot_De        = 1;
flag_plot_ED        = 0;
flag_plot_BD        = 1;

flag_save_fig       = 0;

flag_delete_De      = 1;
flag_delete_ED      = 1;
flag_delete_BD      = 1;

flag_self_weight    = 1;

arrnode  = [2+21+21,2+21+2*21,2+21+3*21]; % blade wing tip nodes

% column for results of wing tip displacement
arrcolumn = [];
nplotDe   = 0; nplotBD = 0; nplotED = 0; 

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
        q = load([model.strSimName '_q.dres']);
        if flag_self_weight == 1
%             copyfile([model.strSimName '_model.m'], ['m' model.strSimName '_model.m']);
%             run(['m' model.strSimName '_model'])
            q0 = load(['self_weight_q.dres']); q = [q0(1,:);q(1:end-1,:)];
        end
        res_de  = fun_extract_wt_tip_deformation(q,time,arrnode,2); 
        save(['res_' var_name '_de'],'res_de');
    end
    nplotDe = length(arrnode);
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
        if strncmp('TipDxb1',strInpVarED{1}(i),11) == 1; arrcolumnED(1,1) = i; end
        if strncmp('TipDyb1',strInpVarED{1}(i),11) == 1; arrcolumnED(1,2) = i; end
        if strncmp('TipDzb1',strInpVarED{1}(i),11) == 1; arrcolumnED(1,3) = i; end
        if strncmp('TipDxb2',strInpVarED{1}(i),11) == 1; arrcolumnED(2,1) = i; end
        if strncmp('TipDyb2',strInpVarED{1}(i),11) == 1; arrcolumnED(2,2) = i; end
        if strncmp('TipDzb2',strInpVarED{1}(i),11) == 1; arrcolumnED(2,3) = i; end
        if strncmp('TipDxb3',strInpVarED{1}(i),11) == 1; arrcolumnED(3,1) = i; end
        if strncmp('TipDyb3',strInpVarED{1}(i),11) == 1; arrcolumnED(3,2) = i; end
        if strncmp('TipDzb3',strInpVarED{1}(i),11) == 1; arrcolumnED(3,3) = i; end
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
        if strncmp('B1TipTDxr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(1,1) = i; end
        if strncmp('B1TipTDyr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(1,2) = i; end
        if strncmp('B1TipTDzr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(1,3) = i; end
        if strncmp('B2TipTDxr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(2,1) = i; end
        if strncmp('B2TipTDyr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(2,2) = i; end
        if strncmp('B2TipTDzr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(2,3) = i; end
        if strncmp('B3TipTDxr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(3,1) = i; end
        if strncmp('B3TipTDyr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(3,2) = i; end
        if strncmp('B3TipTDzr',strInpVarBD{1}(i),11) == 1; arrcolumnBD(3,3) = i; end
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
    var_name = [var_name '_blade' num2str(i)]
    if flag_plot_time_hist == 1
        
%% Plotting time series
        strWhatToComp = []; strleg = {};
        handlesPlot1 = []; handlesPlot2 = []; handlesPlot3 = [];

        figure(); hold on; grid on;
        subplot(3,1,1); hold on; grid on; title(['Tip defl. flap-wise of blade ' num2str(i)]); ylabel('u_1 in m'); xlabel('time in s');
        subplot(3,1,2); hold on; grid on; title(['Tip defl. edge-wise of blade ' num2str(i)]); ylabel('u_2 in m'); xlabel('time in s');
        subplot(3,1,3); hold on; grid on; title(['Tip defl. longit. of blade ' num2str(i)]); ylabel('u_3 in m'); xlabel('time in s');

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];

            subplot(3,1,1); plot(time,res_de(i).u(:,1),'-b','linewidth',lineSize);
            handlePlotDe1 = plot(time(inztime2plotDe),res_de(i).u(inztime2plotDe,1),'.','color','k','marker','o','markersize',markerSize);

            subplot(3,1,2); plot(time,res_de(i).u(:,2),'-b','linewidth',lineSize);
            handlePlotDe2 = plot(time(inztime2plotDe),res_de(i).u(inztime2plotDe,2),'.','color','k','marker','o','markersize',markerSize);

            subplot(3,1,3); plot(time,-res_de(i).u(:,3),'-b','linewidth',lineSize);
            handlePlotDe3 = plot(time(inztime2plotDe),-res_de(i).u(inztime2plotDe,3),'.','color','k','marker','o','markersize',markerSize);

            handlesPlot1 = [handlesPlot1,handlePlotDe1]; 
            handlesPlot2 = [handlesPlot2,handlePlotDe2];
            handlesPlot3 = [handlesPlot3,handlePlotDe3];

            strleg = [strleg,{'DeSiO'}];
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];

            subplot(3,1,1); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,1)),':r','linewidth',lineSize);
            handlePlotED1 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,1)),'.','color','k','marker','+','markersize',markerSize);

            subplot(3,1,2); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,2)),':r','linewidth',lineSize);
            handlePlotED2 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,2)),'.','color','k','marker','+','markersize',markerSize);

            subplot(3,1,3); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,3)),':r','linewidth',lineSize);
            handlePlotED3 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,3)),'.','color','k','marker','+','markersize',markerSize);

            handlesPlot1 = [handlesPlot1,handlePlotED1]; 
            handlesPlot2 = [handlesPlot2,handlePlotED2];
            handlesPlot3 = [handlesPlot3,handlePlotED3];

            strleg = [strleg,{'OpenFAST ED'}];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];

            subplot(3,1,1); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,1)),'--g','linewidth',lineSize);
            handlePlotBD1 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,1)),'.','color','k','marker','s','markersize',markerSize);

            subplot(3,1,2); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,2)),'--g','linewidth',lineSize);
            handlePlotBD2 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,2)),'.','color','k','marker','s','markersize',markerSize);

            subplot(3,1,3); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,3)),'--g','linewidth',lineSize);
            handlePlotBD3 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,3)),'.','color','k','marker','s','markersize',markerSize);

            handlesPlot1 = [handlesPlot1,handlePlotBD1]; 
            handlesPlot2 = [handlesPlot2,handlePlotBD2];
            handlesPlot3 = [handlesPlot3,handlePlotBD3];

            strleg = [strleg,{'OpenFAST BD'}];
        end

        if isempty(strleg) ~= 1
            subplot(3,1,1); leg1 = legend(handlesPlot1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(3,1,2); leg2 = legend(handlesPlot2,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(3,1,3); leg3 = legend(handlesPlot3,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
        end

        if flag_save_fig == 1
            fig1 = gcf; frame_h = get(handle(fig1),'JavaFrame'); set(frame_h,'Maximized',1);
            saveFigure(gcf,['compare_' var_name '_' strWhatToComp '.fig']);
            print(['compare_' var_name '_' strWhatToComp],'-dpng', '-r500');
            set(frame_h,'Maximized',0);
        end
    end
    
    if flag_plot_mean_val == 1
%% plotting mean value
        strleg = {}; strWhatToComp = []; arrColor = [];
        u1_m = []; u2_m = []; u3_m = [];

        figure(); 
        subplot(3,1,1); hold on; grid on; title(['Tip defl. flap-wise of blade ' num2str(i)]); ylabel('mean u_1 in m'); 
        subplot(3,1,2); hold on; grid on; title(['Tip defl. edge-wise of blade ' num2str(i)]); ylabel('mean u_2 in m');
        subplot(3,1,3); hold on; grid on; title(['Tip defl. longit. of blade ' num2str(i)]); ylabel('mean u_3 in m');

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];
            u1_m_de = mean(res_de(i).u(:,1)); 
            u2_m_de = mean(res_de(i).u(:,2)); 
            u3_m_de = -mean(res_de(i).u(:,3));
            u1_m = [u1_m,u1_m_de]; u2_m = [u2_m,u2_m_de]; u3_m = [u3_m,u3_m_de];
            arrColor = [arrColor;[1,0,0]]; strleg = [strleg,{'DeSiO'}];
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];

            u1_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,1)));
            u2_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,2)));
            u3_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,3)));
            u1_m = [u1_m,u1_m_ofastED]; u2_m = [u2_m,u2_m_ofastED]; u3_m = [u3_m,u3_m_ofastED];
            arrColor = [arrColor;[0,1,0]]; strleg = [strleg,{'OpenFAST ED'}];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];
            u1_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,1)));
            u2_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,2)));
            u3_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,3)));
            u1_m = [u1_m,u1_m_ofastBD]; u2_m = [u2_m,u2_m_ofastBD]; u3_m = [u3_m,u3_m_ofastBD];
            arrColor = [arrColor;[0,0,1]]; strleg = [strleg,{'OpenFAST BD'}];
        end

        % bar plot for mean values
        if isempty(u1_m) ~= 1; subplot(3,1,1); bar1 = bar(vertcat([1:length(u1_m)],nan(length(u1_m))),vertcat(u1_m,nan(length(u1_m)))); end
        if isempty(u2_m) ~= 1; subplot(3,1,2); bar2 = bar(vertcat([1:length(u2_m)],nan(length(u2_m))),vertcat(u2_m,nan(length(u2_m)))); end
        if isempty(u3_m) ~= 1; subplot(3,1,3); bar3 = bar(vertcat([1:length(u3_m)],nan(length(u3_m))),vertcat(u3_m,nan(length(u3_m)))); end

        if isempty(strleg) ~= 1
            subplot(3,1,1); leg1 = legend(bar1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(3,1,2); leg2 = legend(bar2,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(3,1,3); leg3 = legend(bar3,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
        end

        if flag_save_fig == 1
            fig2 = gcf; frame_h = get(handle(fig2),'JavaFrame'); set(frame_h,'Maximized',1);
            saveFigure(gcf,['compare_' var_name '_mean_' strWhatToComp '.fig']);
            print(['compare_' var_name '_mean_' strWhatToComp],'-dpng', '-r500');
            set(frame_h,'Maximized',0);
        end
    end

end

return