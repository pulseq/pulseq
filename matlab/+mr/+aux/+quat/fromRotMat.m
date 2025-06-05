function q = fromRotMat(R)
%ROTATE convert normalized quaternion q to rotation matrix

if size(R)~=[3 3]
    error('Only a single rotation matrix can be processed');
end

% fix 'almost zero' or 'almost 1' elements of the rotation matrix
roundedR     = round( R );
boolRound    = abs( roundedR - R ) <= eps; % defines which values should be fixed
R(boolRound) = roundedR(boolRound);
            
if all(R == 0)
    error('Empty (or almost empty) matrix provided in place of a rotation matrix');
end

qs = 0.5 * sqrt( max( 0, R(1,1) + R(2,2) + R(3,3) + 1 ));

if abs(qs) <= eps
    sgn_R23 = 1; sgn_R23(-R(2,3) < 0) = -1;
    sgn_R13 = 1; sgn_R13(-R(1,3) < 0) = -1;
    sgn_R12 = 1; sgn_R12(-R(1,2) < 0) = -1;
    
    q = [0 sqrt( max( 0, -0.5 *( R(2,2) + R(3,3) ))) * sgn_R23 ...
           sqrt( max( 0, -0.5 *( R(1,1) + R(3,3) ))) * sgn_R13 ...
           sqrt( max( 0, -0.5 *( R(1,1) + R(2,2) ))) * sgn_R12 ];
else
    q = [qs 0.25 *( R(3,2) - R(2,3) ) / qs ...
            0.25 *( R(1,3) - R(3,1) ) / qs ...
            0.25 *( R(2,1) - R(1,2) ) / qs ];
end

% normalize the quaternion to account for rounding errors in the rotation matrix, etc
q = mr.aux.quat.normalize(q);
