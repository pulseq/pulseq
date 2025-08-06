function rot = makeRotation( varargin )
% makeRotation Create a rotation extension object
% makeRotation( phi ) - rotation about Z (in plane), angle in rad
% makeRotation( phi, theta ) - rotation on a sphere, angles in rad 
% makeRotation( axis, angle ) - rotation about a given 3D vector by a given angle in rad
% makeRotation( quaternion ) - rotation defined by a unit quaternion
% makeRotation( rot_mat ) - rotation defined by a 3x3 rotation matrix
if nargin<1
    error('makeRotation:invalidArguments','makeRotation - invalid arguments: must supply rotation parameter(s)');
end
switch numel(varargin{1})
    case 1
        phi = varargin{1};
        if nargin<2
            theta = 0.0;
        else
            theta = varargin{2};
        end
        assert( ( phi >= -pi ) && ( phi < 2*pi) , 'makeRotation:invalidTheta',...
            'rotation angle phi (%.2f) is invalid. should be within [-pi,2*pi] radians',phi);
        assert( ( theta >= -pi ) && ( theta <= pi) , 'makeRotation:invalidPhi',...
            'rotation angle theta (%.2f) is invalid. should be within [-pi,pi] radians',theta);
        q1=[cos(theta/2) 0 sin(theta/2) 0]; % y axis
        q2=[cos(phi/2) 0 0 sin(phi/2)]; % z axis
        rot.rotQuaternion=mr.aux.quat.multiply(q2,q1); % ok, looks like the order is right
    case 3
        v=varargin{1};
        v=v./sqrt(sum(v.^2));
        assert(nargin>1);
        phi = varargin{2};
        assert( ( abs(phi) >= 0 ) & ( abs(phi) <= pi) , 'makeRotation:invalidPhi',...
            'rotation angle phi (%.2f) is invalid. should be within [0,pi] radians',phi);
        rot.rotQuaternion=[cos(phi/2) sin(phi/2)*v];
    case 4
        rot.rotQuaternion=mr.aux.quat.normalize(varargin{1});
    case 9
        assert(all(size(varargin{1})==[3 3]));
        rot.rotQuaternion=mr.aux.quat.fromRotMat(varargin{1});
    otherwise
        error('unexpected input to makeRotation');
end

rot.type = 'rot3D';

end

