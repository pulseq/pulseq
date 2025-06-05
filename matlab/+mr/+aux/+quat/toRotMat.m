function r = toRotMat(q)
%ROTATE convert normalized quaternion q to rotation matrix

if length(q(:))~=4
    error('Only a single quaternion expressed as a matlab vector of length 4 can be processed');
end

r = [(1 - 2*(q(3)^2 + q(4)^2)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3)); ...
     2*(q(2)*q(3) + q(1)*q(4)), (1 - 2*(q(2)^2 + q(4)^2)), 2*(q(3)*q(4) - q(1)*q(2)); ...
     2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), (1 - 2*(q(2)^2 + q(3)^2))];

end

