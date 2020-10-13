function [x_after_p] = DOCG_linear(x_before_p,U,m)

x_after_p = zeros(m,1);
for i = 1 : m
    if x_before_p(i) >= 0
        x_after_p(i) = -U;
    else
        x_after_p(i) = +U;
    end
end


