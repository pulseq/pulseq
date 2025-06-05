function r = rotate(q,v)
%ROTATE rotate vector v by quaternion q
%   Rotation ov a vector v by a unit quaternion q can be expressed as qvq',
%   where v is temporarily transformed to a quaternion of a form [0 v]. The
%   function below spells it out explicitly to accelerate calculations 

% A simple quaternion product can be calculated for the vector and scalar
% parts of the quaternion as follows: 
%   vec = s1*v2 + s2*v1 + cross(v1,v2)
%   scalar = s1*s2 - dot(v1,v2)
% We redefine variables as s1=sq=q(0), v1=vq=q(2:4), s2=0, v2=v
% The first part of the product can be writtes as
%   vec1= sq*v + cross(vq,v)
%   scl1= -dot(vq,v)
% The second part of the product is then
%   vec2= -scl1*vq + sq*vec1 + cross(vec1,-vq)
%   scl2= scl1*sq - dot(vec1,-vq)
% but actually scl2 sould be discarded so we don't calculate it

% here is what ChatGPT derives
%r1​=(q1^2​+q2^2​−q3^2​−q4^2​)*v1​+2*(q2​*q3​−q1*​q4​)*v2​+2*(q2*​q4​+q1*​q3​)*v3​
%r2=2*(q2*q3+q1*q4)*v1+(q1^2−q2^2+q3^2−q4^2)*v2+2*(q3*q4−q1*q2)*v3
%r3=2*(q2*q4−q1*q3)*v1+2*(q3*q4+q1*q2)*v2+(q1^2−q2^2−q3^2+q4^2)*v3

r = [(q(:,1)^2 + q(:,2)^2 - q(:,3)^2 - q(:,4)^2)*v(:,1) + 2*(q(:,2)*q(:,3) - q(:,1)*q(:,4))*v(:,2) + 2*(q(:,2)*q(:,4) + q(:,1)*q(:,3))*v(:,3) ...
     2*(q(:,2)*q(:,3) + q(:,1)*q(:,4))*v(:,1) + (q(:,1)^2 - q(:,2)^2 + q(:,3)^2 - q(:,4)^2)*v(:,2) + 2*(q(:,3)*q(:,4) - q(:,1)*q(:,2))*v(:,3) ...
     2*(q(:,2)*q(:,4) - q(:,1)*q(:,3))*v(:,1) + 2*(q(:,3)*q(:,4) + q(:,1)*q(:,2))*v(:,2) + (q(:,1)^2 - q(:,2)^2 - q(:,3)^2 + q(:,4)^2)*v(:,3) ];

end

