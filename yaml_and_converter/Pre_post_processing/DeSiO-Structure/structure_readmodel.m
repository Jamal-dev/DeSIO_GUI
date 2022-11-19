function model = strucure_readmodel()
% =================================================================================================================
model = [];
fid1 = fopen(['simulationinput_structure.txt'],'r');
fid2 = fopen(['simulationinput.txt'],'r');
fid = -1;
    if fid1 ~= -1
        fid = fid1;
    elseif fid2 ~= -1
        fid = fid2;
    end
    if fid >= 0 
        for i = 1:3
            tline = fgets(fid);
        end
        tline = fgetl(fid);
        [str_n,count] = sscanf(tline,'%s%c');
        str = str_n;
        if count>1
            str_n = textscan(str_n,'%s');
            str = cellstr(str_n{1}(1));
            str = str{1};
        end
        model.strSimName = str;
        fclose(fid);
    end

% =================================================================================================================
return