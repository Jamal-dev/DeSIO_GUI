function output_arr = fun_lin_interpolation(input_arr,xhi)
    for i = 1:length(xhi)
        [inz0] = find( xhi(i) <= input_arr(:,1));
        if inz0(1) ~= 1; 
            xhi1 = input_arr(inz0(1)-1,1); y1 = input_arr(inz0(1)-1,2);
            xhi2 = input_arr(inz0(1),1);   y2 = input_arr(inz0(1),2);
            val = (y2-y1)/(xhi2-xhi1)*(xhi(i)-xhi1)+y1;
        else
            val = input_arr(inz0(1),2);
        end
        output_arr(i,1) = [val];
    end
return