import os
import matplotlib.pyplot as plt
import fun_load_file

def uvlm_plot_CL_time():
    os.system('cls||clear')
    # addpath('..\..\files_for_prepostprocessing\DeSiO-Aero\')
    
    currDir = cd
    strfilename = ''
    
    # initializing figures
	f1 = plt.figure()
	plt.grid(visible = True, linewidth = 1)
	plt.xlabel('time s')
	plt.ylabel('C_L')
    # calculate and plot lift coeficient
    model = uvlm_readmodel()
    model.windrotax = [1,0,0]
    t   = fun_load_file(model["strSimName"] + '_uvlm_t.dres')
	t = t[:,1]
    dp  = fun_load_file(model["strSimName"] + '_uvlm_dps_cp.dres')
    qs  = fun_load_file(model["strSimName"] + '_uvlm_qs_nodal.dres')
    CL = uvlm_CL_surface(model,t,qs,dp)
    plt.plot(t,CL,'go--',linewidth = 2, markersize = 20,color = 'b')
    plt.legend('C_L')
    plt.setp(plt.gca.patch, fontweight = 'bold',fontsize = 12)
    # I don't understand the first argument in the below line
	print(strfilename + '_CL_vs_time','-dpng', '-r500')
    