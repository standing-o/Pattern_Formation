function y = lap(s, h)
y = (s(1:end-2, 2:end-1) + s(3:end, 2:end-1) + s(2:end-1, 1:end-2) + s(2:end-1, 3:end) - 4*s(2:end-1, 2:end-1))./h^2;
