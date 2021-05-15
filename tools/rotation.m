function rows = rotation(axis, angle);
s = sin(angle);
c = cos(angle);
if axis == 1
    rows = [1  0 0; 0  c  s;    0 -s  c];
elseif axis == 2
    rows = [c 0 -s; 0  1  0; s -0  c];
else
    rows = [+c  s  0; -s  c  0; +0  0  1];
end
end

