function compare_stress_resultants()
% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
close all
format short

addpath(genpath('U:\PhD-2021-12_Christian_Hente\Pre_post_processing'));

currDir        = cd;
path_desio    = ['\DeSiO\15MW_pitch10_vel0_nspan3\DeSiO'];
path_ofast_BD = ['\OpenFAST_BeamDyn'];
path_ofast_ED = ['\OpenFAST_ElastoDyn'];

str_ofast_BD_filename = 'IEA-15-240-RWT-Onshore';
str_ofast_ED_filename = 'IEA-15-240-RWT-Onshore';

var_name    = 'stress_resultants';

lineSize   = 2;
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

ndeltafactor = 1/3;

arrelem = [50+1,50+1*50+1,50+2*50+1]; % blade wing tip nodes

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
        stress = load([model.strSimName '_stress.dres']);
        res_stress_de  = fun_extract_stress_resulants(stress,time,arrelem); 
        save(['res_' var_name '_de'],'res_stress_de');
    end
    nplotDe        = size(res_stress_de,2);
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
        if strncmp('RootFxb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,1) = i; end
        if strncmp('RootFyb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,2) = i; end
        if strncmp('RootFzb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,3) = i; end
        if strncmp('RootMxb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,4) = i; end
        if strncmp('RootMyb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,5) = i; end
        if strncmp('RootMzb1',strInpVarED{1}(i),8) == 1; arrcolumnED(1,6) = i; end

        if strncmp('RootFxb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,1) = i; end
        if strncmp('RootFyb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,2) = i; end
        if strncmp('RootFzb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,3) = i; end
        if strncmp('RootMxb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,4) = i; end
        if strncmp('RootMyb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,5) = i; end
        if strncmp('RootMzb2',strInpVarED{1}(i),8) == 1; arrcolumnED(2,6) = i; end

        if strncmp('RootFxb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,1) = i; end
        if strncmp('RootFyb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,2) = i; end
        if strncmp('RootFzb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,3) = i; end
        if strncmp('RootMxb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,4) = i; end
        if strncmp('RootMyb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,5) = i; end
        if strncmp('RootMzb3',strInpVarED{1}(i),8) == 1; arrcolumnED(3,6) = i; end
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
        if strncmp('B1N001_FxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,1) = i; end
        if strncmp('B1N001_FyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,2) = i; end
        if strncmp('B1N001_FzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,3) = i; end
        if strncmp('B1N001_MxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,4) = i; end
        if strncmp('B1N001_MyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,5) = i; end
        if strncmp('B1N001_MzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(1,6) = i; end
        
        if strncmp('B2N001_FxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,1) = i; end
        if strncmp('B2N001_FyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,2) = i; end
        if strncmp('B2N001_FzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,3) = i; end
        if strncmp('B2N001_MxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,4) = i; end
        if strncmp('B2N001_MyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,5) = i; end
        if strncmp('B2N001_MzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(2,6) = i; end

        if strncmp('B3N001_FxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,1) = i; end
        if strncmp('B3N001_FyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,2) = i; end
        if strncmp('B3N001_FzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,3) = i; end
        if strncmp('B3N001_MxL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,4) = i; end
        if strncmp('B3N001_MyL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,5) = i; end
        if strncmp('B3N001_MzL',strInpVarBD{1}(i),10) == 1; arrcolumnBD(3,6) = i; end
        
%         if strncmp('B1RootFxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,1) = i; end
%         if strncmp('B1RootFyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,2) = i; end
%         if strncmp('B1RootFzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,3) = i; end
%         if strncmp('B1RootMxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,4) = i; end
%         if strncmp('B1RootMyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,5) = i; end
%         if strncmp('B1RootMzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(1,6) = i; end
% 
%         if strncmp('B2RootFxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,1) = i; end
%         if strncmp('B2RootFyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,2) = i; end
%         if strncmp('B2RootFzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,3) = i; end
%         if strncmp('B2RootMxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,4) = i; end
%         if strncmp('B2RootMyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,5) = i; end
%         if strncmp('B2RootMzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(2,6) = i; end
% 
%         if strncmp('B3RootFxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,1) = i; end
%         if strncmp('B3RootFyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,2) = i; end
%         if strncmp('B3RootFzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,3) = i; end
%         if strncmp('B3RootMxr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,4) = i; end
%         if strncmp('B3RootMyr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,5) = i; end
%         if strncmp('B3RootMzr',strInpVarBD{1}(i),9) == 1; arrcolumnBD(3,6) = i; end
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
        handlesPlot4 = []; handlesPlot5 = []; handlesPlot6 = [];

        figure(); hold on; grid on;
        subplot(2,3,1); hold on; grid on; title(['F_1 of blade ' num2str(i)]); ylabel('F_1 in m'); xlabel('time in s');
        subplot(2,3,2); hold on; grid on; title(['F_2 of blade ' num2str(i)]); ylabel('F_2 in m'); xlabel('time in s');
        subplot(2,3,3); hold on; grid on; title(['F_3 of blade ' num2str(i)]); ylabel('F_3 in m'); xlabel('time in s');
        subplot(2,3,4); hold on; grid on; title(['M_1 of blade ' num2str(i)]); ylabel('M_1 in Nm'); xlabel('time in s');
        subplot(2,3,5); hold on; grid on; title(['M_2 of blade ' num2str(i)]); ylabel('M_2 in Nm'); xlabel('time in s');
        subplot(2,3,6); hold on; grid on; title(['M_3 of blade ' num2str(i)]); ylabel('M_3 in Nm'); xlabel('time in s');

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];

            subplot(2,3,1); plot(time,res_stress_de(i).F1,'-b','linewidth',lineSize);
            handlePlotDe1 = plot(time(inztime2plotDe),res_stress_de(i).F1(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            subplot(2,3,2); plot(time,res_stress_de(i).F2,'-b','linewidth',lineSize);
            handlePlotDe2 = plot(time(inztime2plotDe),res_stress_de(i).F2(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            subplot(2,3,3); plot(time,res_stress_de(i).F3,'-b','linewidth',lineSize);
            handlePlotDe3 = plot(time(inztime2plotDe),res_stress_de(i).F3(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            subplot(2,3,4); plot(time,res_stress_de(i).M1,'-b','linewidth',lineSize);
            handlePlotDe4 = plot(time(inztime2plotDe),res_stress_de(i).M1(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            subplot(2,3,5); plot(time,res_stress_de(i).M2,'-b','linewidth',lineSize);
            handlePlotDe5 = plot(time(inztime2plotDe),res_stress_de(i).M2(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            subplot(2,3,6); plot(time,res_stress_de(i).M3,'-b','linewidth',lineSize);
            handlePlotDe6 = plot(time(inztime2plotDe),res_stress_de(i).M3(inztime2plotDe),'.','color','k','marker','o','markersize',markerSize);

            strleg = [strleg,{'DeSiO'}];

            handlesPlot1 = [handlesPlot1,handlePlotDe1]; handlesPlot2 = [handlesPlot2,handlePlotDe2];
            handlesPlot3 = [handlesPlot3,handlePlotDe3]; handlesPlot4 = [handlesPlot4,handlePlotDe4];
            handlesPlot5 = [handlesPlot5,handlePlotDe5]; handlesPlot6 = [handlesPlot6,handlePlotDe6];        
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];

            subplot(2,3,1); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,1)).*1000,':r','linewidth',lineSize);
            handlePlotED1 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,1)).*1000,'.','color','k','marker','+','markersize',markerSize);

            subplot(2,3,2); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,2)).*1000,':r','linewidth',lineSize);
            handlePlotED2 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,2)).*1000,'.','color','k','marker','+','markersize',markerSize);

            subplot(2,3,3); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,3)).*1000,':r','linewidth',lineSize);
            handlePlotED3 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,3)).*1000,'.','color','k','marker','+','markersize',markerSize);

            subplot(2,3,4); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,4)).*1000,':r','linewidth',lineSize);
            handlePlotED4 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,4)).*1000,'.','color','k','marker','+','markersize',markerSize);

            subplot(2,3,5); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,5)).*1000,':r','linewidth',lineSize);
            handlePlotED5 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,5)).*1000,'.','color','k','marker','+','markersize',markerSize);

            subplot(2,3,6); plot(ofastED_res(:,1),ofastED_res(:,arrcolumnED(i,6)).*1000,':r','linewidth',lineSize);
            handlePlotED6 = plot(ofastED_res(inztime2plotED,1),ofastED_res(inztime2plotED,arrcolumnED(i,6)).*1000,'.','color','k','marker','+','markersize',markerSize);

            strleg = [strleg,{'OpenFAST ED'}];

            handlesPlot1 = [handlesPlot1,handlePlotED1]; handlesPlot2 = [handlesPlot2,handlePlotED2];
            handlesPlot3 = [handlesPlot3,handlePlotED3]; handlesPlot4 = [handlesPlot4,handlePlotED4];
            handlesPlot5 = [handlesPlot5,handlePlotED5]; handlesPlot6 = [handlesPlot6,handlePlotED6];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];

            subplot(2,3,1); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,1)),'--g','linewidth',lineSize);
            handlePlotBD1 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,1)),'.','color','k','marker','s','markersize',markerSize);

            subplot(2,3,2); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,2)),'--g','linewidth',lineSize);
            handlePlotBD2 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,2)),'.','color','k','marker','s','markersize',markerSize);

            subplot(2,3,3); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,3)),'--g','linewidth',lineSize);
            handlePlotBD3 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,3)),'.','color','k','marker','s','markersize',markerSize);

            subplot(2,3,4); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,4)),'--g','linewidth',lineSize);
            handlePlotBD4 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,4)),'.','color','k','marker','s','markersize',markerSize);

            subplot(2,3,5); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,5)),'--g','linewidth',lineSize);
            handlePlotBD5 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,5)),'.','color','k','marker','s','markersize',markerSize);

            subplot(2,3,6); plot(ofastBD_res(:,1),ofastBD_res(:,arrcolumnBD(i,6)),'--g','linewidth',lineSize);
            handlePlotBD6 = plot(ofastBD_res(inztime2plotBD,1),ofastBD_res(inztime2plotBD,arrcolumnBD(i,6)),'.','color','k','marker','s','markersize',markerSize);

            strleg = [strleg,{'OpenFAST BD'}];

            handlesPlot1 = [handlesPlot1,handlePlotBD1]; handlesPlot2 = [handlesPlot2,handlePlotBD2];
            handlesPlot3 = [handlesPlot3,handlePlotBD3]; handlesPlot4 = [handlesPlot4,handlePlotBD4];
            handlesPlot5 = [handlesPlot5,handlePlotBD5]; handlesPlot6 = [handlesPlot6,handlePlotBD6];

        end

        if isempty(strleg) ~= 1
            subplot(2,3,1); leg1 = legend(handlesPlot1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,2); leg2 = legend(handlesPlot2,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,3); leg3 = legend(handlesPlot3,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,4); leg4 = legend(handlesPlot4,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,5); leg5 = legend(handlesPlot5,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,6); leg6 = legend(handlesPlot6,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
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
        F1_m = []; F2_m = []; F3_m = []; M1_m = []; M2_m = []; M3_m = [];

        figure(); 
        subplot(2,3,1); hold on; grid on; title(['F_1 of blade ' num2str(i)]); ylabel('mean F_1 in m'); 
        subplot(2,3,2); hold on; grid on; title(['F_2 of blade ' num2str(i)]); ylabel('mean F_2 in m');
        subplot(2,3,3); hold on; grid on; title(['F_3 of blade ' num2str(i)]); ylabel('mean F_3 in m');
        subplot(2,3,4); hold on; grid on; title(['M_1 of blade ' num2str(i)]); ylabel('mean M_1 in Nm');
        subplot(2,3,5); hold on; grid on; title(['M_2 of blade ' num2str(i)]); ylabel('mean M_2 in Nm');
        subplot(2,3,6); hold on; grid on; title(['M_3 of blade ' num2str(i)]); ylabel('mean M_3 in Nm');

        % DeSiO
        if flag_plot_De == 1
            strWhatToComp = [strWhatToComp '_DeSiO' ];
            F1_m_de = mean(res_stress_de(i).F1); F2_m_de = mean(res_stress_de(i).F2); F3_m_de = mean(res_stress_de(i).F3);
            M1_m_de = mean(res_stress_de(i).M1); M2_m_de = mean(res_stress_de(i).M2); M3_m_de = mean(res_stress_de(i).M3);
            F1_m = [F1_m,F1_m_de]; F2_m = [F2_m,F2_m_de]; F3_m = [F3_m,F3_m_de]; M1_m = [M1_m,M1_m_de]; M2_m = [M2_m,M2_m_de]; M3_m = [M3_m,M3_m_de];
            arrColor = [arrColor;[1,0,0]]; strleg = [strleg,{'DeSiO'}];
        end

        % OpenFAST ElastoDyn
        if flag_plot_ED == 1
            strWhatToComp = [strWhatToComp '_OpenFASTED' ];

            F1_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,1)).*1000);
            F2_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,2)).*1000);
            F3_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,3)).*1000);
            M1_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,4)).*1000);
            M2_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,5)).*1000);
            M3_m_ofastED = mean(ofastED_res(:,arrcolumnED(i,6)).*1000);
            F1_m = [F1_m,F1_m_ofastED]; F2_m = [F2_m,F2_m_ofastED]; F3_m = [F3_m,F3_m_ofastED];
            M1_m = [M1_m,M1_m_ofastED]; M2_m = [M2_m,M2_m_ofastED]; M3_m = [M3_m,M3_m_ofastED];
            arrColor = [arrColor;[0,1,0]]; strleg = [strleg,{'OpenFAST ED'}];
        end

        % OpenFAST BeamDyn
        if flag_plot_BD == 1
            strWhatToComp = [strWhatToComp '_OpenFASTBD' ];

            F1_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,1)));
            F2_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,2)));
            F3_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,3)));
            M1_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,4)));
            M2_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,5)));
            M3_m_ofastBD = mean(ofastBD_res(:,arrcolumnBD(i,6)));
            F1_m = [F1_m,F1_m_ofastBD]; F2_m = [F2_m,F2_m_ofastBD]; F3_m = [F3_m,F3_m_ofastBD];
            M1_m = [M1_m,M1_m_ofastBD]; M2_m = [M2_m,M2_m_ofastBD]; M3_m = [M3_m,M3_m_ofastBD];
            arrColor = [arrColor;[0,0,1]]; strleg = [strleg,{'OpenFAST BD'}];
        end

        % bar plot for mean values
        if isempty(F1_m) ~= 1; subplot(2,3,1); bar1 = bar(vertcat([1:length(F1_m)],nan(length(F1_m))),vertcat(F1_m,nan(length(F1_m)))); end
        if isempty(F2_m) ~= 1; subplot(2,3,2); bar2 = bar(vertcat([1:length(F2_m)],nan(length(F2_m))),vertcat(F2_m,nan(length(F2_m)))); end
        if isempty(F3_m) ~= 1; subplot(2,3,3); bar3 = bar(vertcat([1:length(F3_m)],nan(length(F3_m))),vertcat(F3_m,nan(length(F3_m)))); end
        if isempty(M1_m) ~= 1; subplot(2,3,4); bar4 = bar(vertcat([1:length(M1_m)],nan(length(M1_m))),vertcat(M1_m,nan(length(M1_m)))); end
        if isempty(M2_m) ~= 1; subplot(2,3,5); bar5 = bar(vertcat([1:length(M2_m)],nan(length(M2_m))),vertcat(M2_m,nan(length(M2_m)))); end
        if isempty(M3_m) ~= 1; subplot(2,3,6); bar6 = bar(vertcat([1:length(M3_m)],nan(length(M3_m))),vertcat(M3_m,nan(length(M3_m)))); end

        if isempty(strleg) ~= 1
            subplot(2,3,1); leg1 = legend(bar1,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,2); leg2 = legend(bar2,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,3); leg3 = legend(bar3,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,4); leg4 = legend(bar4,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,5); leg5 = legend(bar5,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
            subplot(2,3,6); leg6 = legend(bar6,strleg); set(gca,'fontsize',FontSize,'fontweight','bold');
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