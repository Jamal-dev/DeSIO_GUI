% function model = fun_read_yaml(str_directory)
function model = fun_read_yaml(str_directory)
    model = [];
    fid = fopen(str_directory,'r');
    if fid == -1
        disp('File .yaml not found!');
        return
    end
    str_obj = []; ilevel = 0; model = []; struct_level = 0;
    while ~feof(fid)
        tline = fgetl(fid);
        if isempty(tline) || ~isempty(strfind(tline,'!!'))
            str_obj      = [];
            ilevel       = 0;
            struct_level = 0;
        else
            bol_stop  = 0;
            ilevel = find_tabs(tline);
            while bol_stop ~= 1
                [bol_stop,str_obj,str_val,ilevel] = fun_check_next_level(fid,tline,bol_stop,str_obj,ilevel);
            end
            
            k_obj = strfind(str_obj,'-');
            str_obj_n = str_obj;
            if ~isempty(k_obj) | struct_level~=0
                if ~isempty(k_obj)
                    [k] = strfind(str_obj,':');
                    struct_level = struct_level + 1;
                    str_obj_o = str_obj(1:k(ilevel));
                    ko = k(ilevel);
                end
                str_obj = [str_obj_o(1:ko-1) '(' num2str(struct_level) ')' str_obj(ko:end)];
                str_obj = strrep(str_obj,'-','');
            end
            eval(['model.' strrep(str_obj,':','.') '=' str_val ';']);
            str_obj = str_obj_n;
        end
    end
return

% recursive function to read model objects and properties
function [bol_stop,str_obj,str_val,ilevel] = fun_check_next_level(fid,str_input,bol_stop,str_obj,ilevel)
    [k] = strfind(str_input,':');
    len_str = length(str_input(k(end)+1:end));
    if len_str~=0
        [g] = strfind(str_obj,':');
        if ilevel == 0
            str_obj = [str_obj strtrim(strjust(str_input(1:k(end)-1),'left'))];
        else
            str_obj = [str_obj(1:g(ilevel)) strtrim(strjust(str_input(1:k(end)-1),'left'))];
        end
        str_val = strjust(str_input(k(end)+1:end),'left');
        bol_stop = 1;
        return
    else
        str_input = [strtrim(str_input) strtrim(strjust(fgetl(fid),'left'))];
        [bol_stop,str_obj,str_val,ilevel] = fun_check_next_level(fid,str_input,bol_stop,str_obj,ilevel);
        ilevel = ilevel + 1;
        if bol_stop == 1
            return
        end
    end
return