function surface = uvlm_readsurfaceinput()
% =================================================================================================================
fid1 = fopen(['surfaceinput.txt'],'r');
fid = -1;
    if fid1 ~= -1
        fid = fid1;
    end
    if fid >= 0 
        for i = 1:3
            tline = fgetl(fid);
        end
        numsurf = str2num(fgetl(fid));
        for k = 1:numsurf
            for i = 1:3
                tline = fgetl(fid);
            end
            tline = fgetl(fid);
            a = sscanf(tline,'%i %i %i %i',4);
            nn = a(1); ne = a(2); nnx = a(3); nny = a(4);
            for i = 1:3
                tline = fgetl(fid);
            end
            for i = 1:nn
                tline = fgetl(fid);
                tline = strrep(tline,'d','e');
                surface(k).coord(i,1:3) = sscanf(tline,'%f %f %f',3);
            end
            for i = 1:3
                tline = fgetl(fid);
            end
            for i = 1:ne
                tline = fgetl(fid);
                surface(k).inz(i,1:4) = sscanf(tline,'%i %i %i %i',4);
            end
        end
    end
% =================================================================================================================
return