function x = pmfRND(f,X)
% pmfRND Ouputs random integers for given PMF
% x = pmfRND(f,X) returns mxn matrix whose values are Zipf distributed,
%     where f > 0 is the PMF (discrete PDF) such that sum(f) = 1;
%     X = [m,n] is the size of desired ouput random matrix.
% if sum(f)>1,error('check input PMF'),end
m = X(1);n = X(2);% size of output random matrix
x = ones(m,n);
u = rand(m,n);
for i = 1:m
    for j = 1:n
        
        p = f(1);
        c = p;%cdf        
        while u(i,j) > c
            x(i,j) = x(i,j) + 1;
            p = f(x(i,j));
            c = c + p;%cdf
        end
    end
end