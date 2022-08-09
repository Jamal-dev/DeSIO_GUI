% =================================================================================================================
% Function to extract geometry from WindIO for structual mesh in DeSiO-Format.
% 
% Author: Christian Hente
% Date: 05.05.2022
%
% input:
%   strName - type of surface of (airfoil) cross-section
%   model - component beam object specified in WindIO
%   materials - material data specified in WindIO
%   scale_opt - scaling factor for scaling in longitudinal direction (optional input)
% output:
%     beam - beam object containing coordinates and connectivity for creating structural mesh
% =================================================================================================================
function beam = fun_extract_beam_data(strName,model_beam,materials)
% =================================================================================================================
    % span-wise discretization
    beam.M = model_beam.M_struc;

    % initializing stiffness and mass matrix
    beam.arr_stiff_matrix = zeros(beam.M,21);
    beam.arr_mass_matrix = zeros(beam.M,6);
    O = 0*ones(beam.M,1);

    % nodal natural coordinates in span-wise direction
    xhi_x  = [0:1/beam.M:1];
    
    % element natural coordinates in span-wise direction
    xhi_x_elem = (xhi_x(1:end-1)+xhi_x(2:end))/2;

    % interpolating coordinates of reference axis according to span-wise discretization
    beam.arr_xre_x = fun_lin_interpolation([model_beam.reference_axis.x.grid',model_beam.reference_axis.x.values'],xhi_x);
    beam.arr_xre_y = fun_lin_interpolation([model_beam.reference_axis.y.grid',model_beam.reference_axis.y.values'],xhi_x);
    beam.arr_xre_z = fun_lin_interpolation([model_beam.reference_axis.z.grid',model_beam.reference_axis.z.values'],xhi_x);
    
    % interpolating twist angle according to span-wise discretization
    if isfield(model_beam,'twist');
        beam.arr_twist = fun_lin_interpolation([model_beam.twist.grid',model_beam.twist.values'],xhi_x);
    end
    
    % if-statement according to type of surface cross-section
    if strfind(strName,'blade') ~= 0
        % interpolating along span-wise direction
        if isfield(model_beam.elastic_properties_mb,'six_x_six')
            beam.arr_stiff_matrix = [];
            beam.arr_mass_matrix  = [];
            % stiffness and mass/interia terms should be already converted to DeSiO-Format, i.e. Voigt notation
            for i = 1:21
                beam.arr_stiff_matrix(:,i) = fun_lin_interpolation([model_beam.elastic_properties_mb.six_x_six.stiff_matrix.grid',model_beam.elastic_properties_mb.six_x_six.stiff_matrix.values(:,i)],xhi_x_elem);
            end
            for i = 1:6
                beam.arr_mass_matrix(:,i)  = fun_lin_interpolation([model_beam.elastic_properties_mb.six_x_six.inertia_matrix.grid',model_beam.elastic_properties_mb.six_x_six.inertia_matrix.values(:,i)],xhi_x_elem);
            end
            beam.dissipation(1:2) = [model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
        else
            
            beam.arr_EA  = zeros(beam.M,1); beam.arr_GA1 = zeros(beam.M,1); beam.arr_GA2 = zeros(beam.M,1);
            beam.arr_EI1 = zeros(beam.M,1); beam.arr_EI2 = zeros(beam.M,1); beam.arr_GI3 = zeros(beam.M,1);
            beam.arr_ES1 = zeros(beam.M,1); beam.arr_ES2 = zeros(beam.M,1); beam.arr_EI12 = zeros(beam.M,1);
            beam.arr_GS1 = zeros(beam.M,1); beam.arr_GS2 = zeros(beam.M,1);
            
            beam.arr_rhoA  = zeros(beam.M,1); beam.arr_rhoI1 = zeros(beam.M,1); beam.arr_rhoI2  = zeros(beam.M,1);
            beam.arr_rhoS1 = zeros(beam.M,1); beam.arr_rhoS2 = zeros(beam.M,1); beam.arr_rhoI12 = zeros(beam.M,1);
            
            if isfield(model_beam.elastic_properties_mb,'EA')
                beam.arr_EA   = fun_lin_interpolation([model_beam.elastic_properties_mb.EA.grid',model_beam.elastic_properties_mb.EA.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'GA1')
                beam.arr_GA1  = fun_lin_interpolation([model_beam.elastic_properties_mb.GA1.grid',model_beam.elastic_properties_mb.GA1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'GA2')
                beam.arr_GA2  = fun_lin_interpolation([model_beam.elastic_properties_mb.GA2.grid',model_beam.elastic_properties_mb.GA2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI1')
                beam.arr_EI1  = fun_lin_interpolation([model_beam.elastic_properties_mb.EI1.grid',model_beam.elastic_properties_mb.EI1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI2')
                beam.arr_EI2  = fun_lin_interpolation([model_beam.elastic_properties_mb.EI2.grid',model_beam.elastic_properties_mb.EI2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'GI3')
                beam.arr_GI3  = fun_lin_interpolation([model_beam.elastic_properties_mb.GI3.grid',model_beam.elastic_properties_mb.GI3.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'ES1')
                beam.arr_ES1  = fun_lin_interpolation([model_beam.elastic_properties_mb.ES1.grid',model_beam.elastic_properties_mb.ES1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'ES2')
                beam.arr_ES2  = fun_lin_interpolation([model_beam.elastic_properties_mb.ES2.grid',model_beam.elastic_properties_mb.ES2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'GS1')
                beam.arr_GS1  = fun_lin_interpolation([model_beam.elastic_properties_mb.GS1.grid',model_beam.elastic_properties_mb.GS1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'GS2')
                beam.arr_GS2  = fun_lin_interpolation([model_beam.elastic_properties_mb.GS2.grid',model_beam.elastic_properties_mb.GS2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI12')
                beam.arr_EI12 = fun_lin_interpolation([model_beam.elastic_properties_mb.EI12.grid',model_beam.elastic_properties_mb.EI12.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoA')
                beam.arr_rhoA   = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoA.grid',model_beam.elastic_properties_mb.rhoA.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI1')
                beam.arr_rhoI1  = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoI1.grid',model_beam.elastic_properties_mb.rhoI1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI2')
            beam.arr_rhoI2  = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoI2.grid',model_beam.elastic_properties_mb.rhoI2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoS1')
                beam.arr_rhoS1  = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoS1.grid',model_beam.elastic_properties_mb.rhoS1.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoS2')
                beam.arr_rhoS2  = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoS2.grid',model_beam.elastic_properties_mb.rhoS2.values'],xhi_x_elem);
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI12')
                beam.arr_rhoI12 = fun_lin_interpolation([model_beam.elastic_properties_mb.rhoI12.grid',model_beam.elastic_properties_mb.rhoI12.values'],xhi_x_elem);
            end
            % GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
            beam.arr_stiff_matrix(:,1:21) = [beam.arr_GA1, beam.arr_GA2, beam.arr_EA, beam.arr_EI1, beam.arr_EI2, beam.arr_GI3, ...
                                            O, O, O, beam.arr_GS2, -beam.arr_GS1, O, O, O, O, O, beam.arr_ES1, -beam.arr_EI12, -beam.arr_ES2, O, O];
            
            % rhoA, rhoI2, rhoI1, rhoI12, rhoS1, rhoS2
            beam.arr_mass_matrix(:,1:6)   = [beam.arr_rhoA, beam.arr_rhoI2, beam.arr_rhoI1, beam.arr_rhoI12, beam.arr_rhoS1, beam.arr_rhoS2];
            beam.dissipation(1:2)         = [model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
        end
    end
    if strfind(strName,'pipe')~=0
        % interpolating cross-section properties along span-wise direction
        E = 0.0; G = 0.0; rho = 0.0;
        if isfield(model_beam,'material')
            strmaterial = model_beam.material;
            for i = 1:size(materials,1)
                if strfind(materials{i}.name,strmaterial)
                    if length(materials{i}.name) == length(strmaterial)
                        E   = materials{i}.E;
                        nu  = materials{i}.nu;
                        rho = materials{i}.rho;
                        G   = E/(2.0*(1.0+nu));
                    end
                end
            end
        end
        if isfield(model_beam,'elastic_properties_mb')
            E   = model_beam.elastic_properties_mb.E;
            G   = E/(2.0*(1+model_beam.elastic_properties_mb.nu));
            rho = model_beam.elastic_properties_mb.rho;
        end

        k1  = 1.0;
        k2  = 1.0;
        if isfield(model_beam,'shear_factor')
            k1 = model_beam.shear_factor.k1;
            k2 = model_beam.shear_factor.k2;
        end

        beam.arr_outer_diameter = fun_lin_interpolation([model_beam.outer_diameter.grid',model_beam.outer_diameter.values'],xhi_x_elem);
        beam.arr_thickness      = fun_lin_interpolation([model_beam.thickness.grid',model_beam.thickness.values'],xhi_x_elem);
        beam.arr_EA             = E*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_GA1            = k1*G*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_GA2            = k2*G*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_EI1            = E*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_EI2            = E*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_GI3            = G*pi/32*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_ES1            = beam.arr_EA*0;
        beam.arr_ES2            = beam.arr_EA*0;
        beam.arr_GS2            = beam.arr_EA*0;
        beam.arr_GS1            = beam.arr_EA*0;
        beam.arr_EI12           = beam.arr_EA*0;

        beam.arr_rhoA           = rho*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_rhoI1          = rho*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_rhoI2          = rho*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_rhoI12         = beam.arr_rhoA*0;
        beam.arr_rhoS1          = beam.arr_rhoA*0;
        beam.arr_rhoS2          = beam.arr_rhoA*0;

        % GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
        beam.arr_stiff_matrix(:,1:21) = [beam.arr_GA1, beam.arr_GA2, beam.arr_EA, beam.arr_EI1, beam.arr_EI2, beam.arr_GI3, ...
                                        O, O, O, beam.arr_GS2, -beam.arr_GS1, O, O, O, O, O, beam.arr_ES1, -beam.arr_EI12, -beam.arr_ES2, O, O];

        beam.arr_mass_matrix(:,1:6) = [beam.arr_rhoA, beam.arr_rhoI2, beam.arr_rhoI1, beam.arr_rhoI12, beam.arr_rhoS1, beam.arr_rhoS2];
        beam.dissipation(1:2)       = [ model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
    end
    
return
