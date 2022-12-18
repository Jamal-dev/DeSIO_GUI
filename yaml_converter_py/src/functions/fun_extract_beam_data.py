from classes import Beam
def fun_extract_beam_data(strName,model_beam,materials):
    """
         input:
           strName - type of surface of (airfoil) cross-section
           model - component beam object specified in WindIO
           materials - material data specified in WindIO
           scale_opt - scaling factor for scaling in longitudinal direction (optional input)
         output:
             beam - beam object containing coordinates and connectivity for creating structural mesh
    """
  
    # span-wise discretization
    beam = Beam(strName,model_beam,materials)
    return beam