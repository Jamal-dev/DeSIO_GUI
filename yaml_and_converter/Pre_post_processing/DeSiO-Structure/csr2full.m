% Function to convert csr (compressed sparse row) matrix format to full (dense) matrix format
function A = csr2full(dA,rA,jA,m,n)

A = zeros(m,n);                                  
for i = 1:m                                      
    for k = rA(i):rA(i+1) -1                     
        j = jA(k);                               
        A(i,j) = dA(k);                          
    end
end


% dA : Value
% rA : Row Pointer
% jA : Column-index


