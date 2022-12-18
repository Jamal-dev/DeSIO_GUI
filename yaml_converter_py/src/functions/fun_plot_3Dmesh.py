import matplotlib.pyplot as plt
import numpy as np
def fun_plot_3Dmesh(str_b,mesh):
    # plot system
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.tight_layout()
    for i  in range(len(mesh)):
        xx = []; yy = []; zz = []
        for i_air in range(mesh[i].connectivity.shape[0]):
            if str_b == '2d' or str_b == '2D':
                x = mesh[i].nodes[mesh[i].connectivity[i_air,:],0]
                y = mesh[i].nodes[mesh[i].connectivity[i_air,:],1]
                z = mesh[i].nodes[mesh[i].connectivity[i_air,:],2]
                # x = x.reshape(2,2)
                # y = y.reshape(2,2)
                # z = z.reshape(2,2)
                xx.append(x)
                yy.append(y)
                zz.append(z)
                # ax.scatter(x,y,z)
                # ax.plot_trisurf(x, y, z,alpha=0.5)
        ax.plot_surface(np.asarray(xx), np.asarray(yy), 
                        np.asarray(zz),alpha=0.5)