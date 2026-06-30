function m = x(al)
% X Rotation matrix for rotation about X to an angle al given in radians
c = cos(al);
s = sin(al);
m = [1 0 0; 0 c -s; 0 s c];
end

