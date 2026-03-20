function tests = testRotate
    tests = functiontests(localfunctions);
end

%% --- mr.rotate standalone tests ---

%% Test identity rotation (angle=0) leaves gradient unchanged
function test_zero_angle(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('z', 0, gx);
    testCase.verifyEqual(length(out), 1);
    testCase.verifyEqual(out{1}.channel, 'x');
    testCase.verifyEqual(out{1}.area, gx.area, 'AbsTol', 1);
end

%% Test 90-degree rotation about z: x -> y
function test_90deg_z_rotation(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('z', pi/2, gx);
    areas = getAreaStruct(out);
    testCase.verifyEqual(abs(areas.y), abs(gx.area), 'AbsTol', 10);
    testCase.verifyLessThan(abs(areas.x), 10, 'X component should vanish');
end

%% Test 90-degree rotation about x: y -> z
function test_90deg_x_rotation(testCase)
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('x', pi/2, gy);
    areas = getAreaStruct(out);
    testCase.verifyEqual(abs(areas.z), abs(gy.area), 'AbsTol', 10);
    testCase.verifyLessThan(abs(areas.y), 10, 'Y component should vanish');
end

%% Test 90-degree rotation about y: z -> x
function test_90deg_y_rotation(testCase)
    gz = mr.makeTrapezoid('z', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('y', pi/2, gz);
    areas = getAreaStruct(out);
    testCase.verifyEqual(abs(areas.x), abs(gz.area), 'AbsTol', 10);
    testCase.verifyLessThan(abs(areas.z), 10, 'Z component should vanish');
end

%% Test gradient parallel to rotation axis is unchanged
function test_parallel_gradient_unchanged(testCase)
    gz = mr.makeTrapezoid('z', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('z', pi/3, gz);
    testCase.verifyEqual(length(out), 1);
    testCase.verifyEqual(out{1}.channel, 'z');
    testCase.verifyEqual(out{1}.area, gz.area, 'AbsTol', 1);
end

%% Test 45-degree rotation splits into two components
function test_45deg_split(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/4;
    out = mr.rotate('z', angle, gx);
    areas = getAreaStruct(out);
    testCase.verifyEqual(areas.x, gx.area * cos(angle), 'AbsTol', 10);
    testCase.verifyEqual(areas.y, gx.area * sin(angle), 'AbsTol', 10);
end

%% Test non-gradient events pass through
function test_non_grad_passthrough(testCase)
    adc = mr.makeAdc(128, 'Duration', 1e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('z', pi/4, gx, adc);
    has_adc = false;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'type') && strcmp(out{i}.type, 'adc')
            has_adc = true;
        end
    end
    testCase.verifyTrue(has_adc, 'ADC event should pass through');
end

%% Test with system parameter
function test_with_system(testCase)
    sys = mr.opts();
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    out = mr.rotate('z', pi/4, 'system', sys, gx);
    testCase.verifyTrue(~isempty(out));
    areas = getAreaStruct(out);
    testCase.verifyEqual(areas.x, gx.area * cos(pi/4), 'AbsTol', 10);
end

%% Test invalid axis throws error
function test_invalid_axis_error(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    verifyErrorThrown(testCase, @() mr.rotate('w', pi/4, gx));
end

%% Test non-scalar angle throws error
function test_non_scalar_angle_error(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    verifyErrorThrown(testCase, @() mr.rotate('z', [pi/4 pi/2], gx));
end

%% Test two gradients on perpendicular axes
function test_two_gradients(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', 'Area', 2000, 'Duration', 2e-3);
    angle = pi/6;
    out = mr.rotate('z', angle, gx, gy);
    areas = getAreaStruct(out);
    % x_out = gx*cos - gy*sin, y_out = gx*sin + gy*cos
    testCase.verifyEqual(areas.x, gx.area*cos(angle) - gy.area*sin(angle), 'AbsTol', 15);
    testCase.verifyEqual(areas.y, gx.area*sin(angle) + gy.area*cos(angle), 'AbsTol', 15);
end

%% --- Cross-validation: mr.rotate vs mr.rotate3D ---

%% Test mr.rotate('z', angle) matches mr.rotate3D(rotz(angle))
function test_cross_validate_z(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', 'Area', 500, 'Duration', 2e-3);
    angle = pi/5;
    out_rotate = mr.rotate('z', angle, gx, gy);
    out_rotate3D = mr.rotate3D(rotz(angle), gx, gy);
    map_r  = getAreaStruct(out_rotate);
    map_3d = getAreaStruct(out_rotate3D);
    testCase.verifyEqual(map_r.x, map_3d.x, 'AbsTol', 5, ...
        'X area should match between rotate and rotate3D for z-axis');
    testCase.verifyEqual(map_r.y, map_3d.y, 'AbsTol', 5, ...
        'Y area should match between rotate and rotate3D for z-axis');
    testCase.verifyEqual(map_r.z, map_3d.z, 'AbsTol', 5, ...
        'Z area should match between rotate and rotate3D for z-axis');
end

%% Test mr.rotate('x', angle) matches mr.rotate3D(rotx(angle))
function test_cross_validate_x(testCase)
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', 'Area', 500, 'Duration', 2e-3);
    angle = pi/7;
    out_rotate = mr.rotate('x', angle, gy, gz);
    out_rotate3D = mr.rotate3D(rotx(angle), gy, gz);
    map_r  = getAreaStruct(out_rotate);
    map_3d = getAreaStruct(out_rotate3D);
    testCase.verifyEqual(map_r.x, map_3d.x, 'AbsTol', 5, ...
        'X area should match between rotate and rotate3D for x-axis');
    testCase.verifyEqual(map_r.y, map_3d.y, 'AbsTol', 5, ...
        'Y area should match between rotate and rotate3D for x-axis');
    testCase.verifyEqual(map_r.z, map_3d.z, 'AbsTol', 5, ...
        'Z area should match between rotate and rotate3D for x-axis');
end

%% Test mr.rotate('y', angle) matches mr.rotate3D(roty(angle))
function test_cross_validate_y(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', 'Area', 500, 'Duration', 2e-3);
    angle = pi/3;
    out_rotate = mr.rotate('y', angle, gx, gz);
    out_rotate3D = mr.rotate3D(roty(angle), gx, gz);
    map_r  = getAreaStruct(out_rotate);
    map_3d = getAreaStruct(out_rotate3D);
    testCase.verifyEqual(map_r.x, map_3d.x, 'AbsTol', 5, ...
        'X area should match between rotate and rotate3D for y-axis');
    testCase.verifyEqual(map_r.y, map_3d.y, 'AbsTol', 5, ...
        'Y area should match between rotate and rotate3D for y-axis');
    testCase.verifyEqual(map_r.z, map_3d.z, 'AbsTol', 5, ...
        'Z area should match between rotate and rotate3D for y-axis');
end

%% Test cross-validation with single gradient on all three axes
function test_cross_validate_single_grad_all_axes(testCase)
    angles = [pi/6, pi/4, pi/3];
    axes_list = {'x', 'y', 'z'};
    rot_funcs = {@rotx, @roty, @rotz};
    for a = 1:3
        for ch = 1:3
            g = mr.makeTrapezoid(axes_list{ch}, 'Area', 1000, 'Duration', 2e-3);
            angle = angles(a);
            out_r  = mr.rotate(axes_list{a}, angle, g);
            out_3d = mr.rotate3D(rot_funcs{a}(angle), g);
            map_r  = getAreaStruct(out_r);
            map_3d = getAreaStruct(out_3d);
            testCase.verifyEqual(map_r.x, map_3d.x, 'AbsTol', 10, ...
                sprintf('X mismatch: rot %s, grad %s', axes_list{a}, axes_list{ch}));
            testCase.verifyEqual(map_r.y, map_3d.y, 'AbsTol', 10, ...
                sprintf('Y mismatch: rot %s, grad %s', axes_list{a}, axes_list{ch}));
            testCase.verifyEqual(map_r.z, map_3d.z, 'AbsTol', 10, ...
                sprintf('Z mismatch: rot %s, grad %s', axes_list{a}, axes_list{ch}));
        end
    end
end

%% --- Helpers ---

function as = getAreaStruct(in)
    as.x = 0; as.y = 0; as.z = 0;
    for i = 1:length(in)
        if isstruct(in{i}) && isfield(in{i}, 'channel')
            as.(in{i}.channel) = in{i}.area;
        end
    end
end

function verifyErrorThrown(testCase, funcHandle)
    didError = false;
    try
        funcHandle();
    catch
        didError = true;
    end
    testCase.verifyTrue(didError, 'Expected an error, but none was thrown.');
end

%% Local rotation matrix constructors (equivalent to rotx/roty/rotz)
function R = rotx(angle)
    R = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
end

function R = roty(angle)
    R = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
end

function R = rotz(angle)
    R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
end
