function m = y(al)
% Y Rotation matrix for rotation about Y to an angle al given in radians
c = cos(al);
s = sin(al);
m = [c -s 0; s c 0; 0 0 1];
end

