function q = conjugate(q)
%CONJUGATE Calculate the conjugate of a quaternion
%  A single quaternion is represented as a 1 x 4 vector with the first
%  component being the real part and the 2nd to 4th components
%  corresponding to the complex vector part.Collections of N quaternions
%  can be stored as N x 4 matrices. 

q(:,2:4) = -q(:,2:4);
end

