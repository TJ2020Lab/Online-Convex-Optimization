%% compute the sub-differential related to the inequalities


function [r] = iclr_sub_ineqs(m,x,L,U);
r = zeros(2*m,1);

for i = 1 : m
    if x(i) < L
        r(i) = -1;
        r(m+i) = 0;
    elseif x(i) > U
        r(i) = 0;
        r(m+i) = 1;
    else
        r(i) = 0;
        r(m+i) = 0;
    end
end

