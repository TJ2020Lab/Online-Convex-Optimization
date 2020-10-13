%% project onto the box


function [r] = box_proj(m,x,L,U);
r = zeros(m,1);

for i = 1 : m
    if x(i) < L
        r(i) = L;
    elseif x(i) > U;
        r(i) = U;
    else
        r(i) = x(i);
    end
end