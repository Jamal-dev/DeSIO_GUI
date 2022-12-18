from classes import UVLM
def fun_extract_uvlm_data(strName,model_uvlm,airfoils):
    """   
     input:
       strName - type of surface of (airfoil) cross-section
       model - component uvlm object specified in WindIO
       airfoils - airfoil data specified in WindIO
     output:
         uvlm_ob - uvlm object containing coordinates and connectivity for creating aerodynamic grid
     =================================================================================================================
    """
    uvlm_ob = UVLM(strName,model_uvlm,airfoils)
    # for debugging
#     for var in vars(uvlm_ob):
#         print(f"{var} : ",getattr(uvlm_ob, var))
#         print('*'*80)
#         print('*'*80)
#         print('*'*80)
    return uvlm_ob