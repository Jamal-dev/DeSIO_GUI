import numpy as np

def uvlm_m_surface(model,t,qs,dp,x0):
    M = []
    for k in range (1,length(t)+1):
        a = 0
        mcp = [0,0,0]
        for i in range (1, model['nsurfaces']+1):
            nodes = model['surfaces'][i]['nodes']
            for j in range (1, model['surfaces'][i]['nelements'] + 1):            
                c   = model['surfaces'][i]['connectivity'][j:,:]                
                n1  = c[1], n2 = c[2], n3 = c[3], n4 = c[4]
                # Area of surface
                # n1,n2... should be indexed or is it function parameter ?!?!?
                diag1 = qs(k,nodes[n1]['indices_q'])-qs(k,nodes[n3]['indices_q'])
                diag2 = qs(k,nodes[n2]['indices_q'])-qs(k,nodes[n4]['indices_q'])
                e1 = (diag2-diag1)
                e2 =-(diag1+diag2)
                e1 = e1/np.norm(e1)
                e2 = e2/np.norm(e2)
                nj = np.cross(e1, e2)
                area1 = 0.5*np.norm(np.cross(qs(k,nodes[n2]['indices_q'])-qs(k,nodes[n1]['indices_q']), qs(k,nodes[n4]['indices_q'])-qs(k,nodes[n1]['indices_q'])))
                area2 = 0.5*np.norm(np.cross(qs(k,nodes[n4]['indices_q'])-qs(k,nodes[n3]['indices_q']), qs(k,nodes[n2]['indices_q'])-qs(k,nodes[n3]['indices_q'])))
                Aj = area1 + area2
                # coordinates of control point
                xcpj = 0.25*np.transpose( qs(k,nodes[n1]['indices_q']) + qs(k,nodes[n2]['indices_q']) + qs(k,nodes[n3]['indices_q']) + qs(k,nodes[n4]['indices_q']) )
                xaj  = xcpj - x0
                # resultant force/dynamic pressure due to pressure
                Fj  = dp(k,j+a)*Aj*np.transpose(nj)
                # mj  = np.cross(Fj,xaj)';
                mj  = np.transpose(np.cross(xaj,Fj))
                mcp = mcp + mj
            # global lift = projecting Cp_vector in lift direction
            a = a + model['surfaces'][i]['nelements']
        M.append(mcp)
    return M