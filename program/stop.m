function [value, isterminal, direction] = stop(t, y)
    value = y(1);
    isterminal = 1;
    direction = 0;
end

