function qout = multiply(q1,q2)
%MULTIPLY Calculate the product of two quaternions.
%  A single quaternion is represented as a 1 x 4 vector with the first
%  component being the real part and the 2nd to 4th components
%  corresponding to the complex vector part.Collections of N quaternions
%  can be stored as N x 4 matrices. 

% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)
vec = [q1(:,1).*q2(:,2) q1(:,1).*q2(:,3) q1(:,1).*q2(:,4)] + ...
         [q2(:,1).*q1(:,2) q2(:,1).*q1(:,3) q2(:,1).*q1(:,4)]+...
         [ q1(:,3).*q2(:,4)-q1(:,4).*q2(:,3) ...
           q1(:,4).*q2(:,2)-q1(:,2).*q2(:,4) ...
           q1(:,2).*q2(:,3)-q1(:,3).*q2(:,2)];

% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q1(:,1).*q2(:,1) - q1(:,2).*q2(:,2) - ...
             q1(:,3).*q2(:,3) - q1(:,4).*q2(:,4);

qout = [scalar  vec];

end

