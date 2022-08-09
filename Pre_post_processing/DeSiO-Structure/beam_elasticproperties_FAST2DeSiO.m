function beam_elasticproperties_FAST2DeSiO()
% =================================================================================================================
clc; clear all; close all; format long;
addpath('G:\10_DeSiO\yaml\yaml\yaml');

% =================================================================================================================
strfilename  = 'IEA-15'
str_dir_yaml = 'IEA-15-240-RWT_test.yaml';    
% =================================================================================================================

% =================================================================================================================
%                   READING WINDIO-FORMAT with informations in blade elastic properties
% =================================================================================================================
    model = YAML.read(str_dir_yaml);
% =================================================================================================================
% reformulate to DeSiO-notation that used Voigt notation

% stiffness matrix
inzstiff  = [1,7,12,16,19,21,20,18,15,11,6,5,4,3,2,8,13,17,14,10,9];
T_row_voigt = zeros(21,21);
for i = 1:length(inzstiff); T_row_voigt(i,inzstiff(i)) = 1; end
arr_stiffness = model.components.blade.elastic_properties_mb.six_x_six.stiff_matrix.values;
for i = 1:size(arr_stiffness,1); arr_stiff_matrix_DeSiO(i,:) = (T_row_voigt*arr_stiffness(i,:)'); end
% mass matrix[6]
inzmass = [1,19,16,-17,13,11];
T_row_voigt = zeros(6,21);
for i = 1:length(inzmass); T_row_voigt(i,abs(inzmass(i))) = 1*sign(inzmass(i)); end
arr_mass = model.components.blade.elastic_properties_mb.six_x_six.inertia_matrix.values;
for i = 1:size(arr_mass,1); arr_mass_matrix_DeSiO(i,:) = (T_row_voigt*arr_mass(i,:)'); end

% print matrices to files
fid = fopen([strfilename '_beamprop_FAST2DeSiO.txt'],'w');
fprintf(fid,'stiffness matrix\n');
for i = 1:size(arr_stiff_matrix_DeSiO,1)
    fprintf(fid,'-  [');
    fprintf(fid,['%10.3e,'],arr_stiff_matrix_DeSiO(i,:));
    fprintf(fid,']\n');
end

fprintf(fid,'\n');
fprintf(fid,'mass matrix\n');
for i = 1:size(arr_mass_matrix_DeSiO,1)
    fprintf(fid,'-  [');
    fprintf(fid,['%10.3e,'],arr_mass_matrix_DeSiO(i,:));
    fprintf(fid,']\n');
end
fclose(fid);

return