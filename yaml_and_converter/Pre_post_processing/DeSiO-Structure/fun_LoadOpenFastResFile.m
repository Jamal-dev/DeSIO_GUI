function [ofast_res,strInpVar] = fun_LoadOpenFastResFile(strfilename)
    q   = [];
    fid = fopen(strfilename);
    ofast_res = cell2mat(textscan(fid,'','headerlines',8));
    fclose(fid);

    fid = fopen(strfilename);
    for i = 1:6; strval = fgetl(fid); end
    strval = fgetl(fid);
    strInpVar = textscan(strval,'%s'); 
    fclose(fid);