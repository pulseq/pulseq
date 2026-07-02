function m = y(al)
% Y Rotation matrix for rotation about Y to an angle al given in radians
c = cos(al);
s = sin(al);
m = [c 0 s; 0 1 0; -s 0 c];
end

