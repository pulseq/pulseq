function tests = testRotate3D
    tests = functiontests(localfunctions);
end

%% Test identity rotation leaves gradients unchanged
function test_identity_rotation(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', 'Area', 2000, 'Duration', 2e-3);
    R = eye(3);
    out = mr.rotate3D(R, gx, gy);
    % Should have 2 gradient events
    testCase.verifyEqual(length(out), 2);
    % Areas should be preserved
    areas = zeros(1, 2);
    for i = 1:2
        areas(i) = out{i}.area;
    end
    testCase.verifyTrue(any(abs(areas - 1000) < 10), 'X area should be preserved');
    testCase.verifyTrue(any(abs(areas - 2000) < 10), 'Y area should be preserved');
end

%% Test 90-degree rotation about z: x->y, y->-x
function test_90deg_z_rotation(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/2;
    R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
    out = mr.rotate3D(R, gx);
    % gx rotated 90° about z should become gy
    testCase.verifyTrue(length(out) >= 1);
    % Find the y-channel gradient
    found_y = false;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel') && out{i}.channel == 'y'
            found_y = true;
            testCase.verifyEqual(abs(out{i}.area), abs(gx.area), 'AbsTol', 10);
        end
    end
    testCase.verifyTrue(found_y, 'Should have a y-channel gradient after 90° z rotation');
end

%% Test 90-degree rotation about x: y->z, z->-y
function test_90deg_x_rotation(testCase)
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/2;
    R = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
    out = mr.rotate3D(R, gy);
    testCase.verifyTrue(length(out) >= 1);
    % gy rotated 90° about x should become gz
    found_z = false;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel') && out{i}.channel == 'z'
            found_z = true;
            testCase.verifyEqual(abs(out{i}.area), abs(gy.area), 'AbsTol', 10);
        end
    end
    testCase.verifyTrue(found_z, 'Should have a z-channel gradient after 90° x rotation');
end

%% Test 90-degree rotation about y: z->x, x->-z
function test_90deg_y_rotation(testCase)
    gz = mr.makeTrapezoid('z', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/2;
    R = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
    out = mr.rotate3D(R, gz);
    testCase.verifyTrue(length(out) >= 1);
    % gz rotated 90° about y should become gx
    found_x = false;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel') && out{i}.channel == 'x'
            found_x = true;
            testCase.verifyEqual(abs(out{i}.area), abs(gz.area), 'AbsTol', 10);
        end
    end
    testCase.verifyTrue(found_x, 'Should have an x-channel gradient after 90° y rotation');
end

%% Test 45-degree rotation about z splits x into x and y
function test_45deg_z_rotation(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/4;
    R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
    out = mr.rotate3D(R, gx);
    % Should produce both x and y components with equal magnitude
    area_x = 0; area_y = 0;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel')
            if out{i}.channel == 'x'
                area_x = out{i}.area;
            elseif out{i}.channel == 'y'
                area_y = out{i}.area;
            end
        end
    end
    expected = gx.area * cos(angle);
    testCase.verifyEqual(area_x, expected, 'AbsTol', 10);
    testCase.verifyEqual(area_y, gx.area * sin(angle), 'AbsTol', 10);
end

%% Test 45-degree rotation about x splits y into y and z
function test_45deg_x_rotation(testCase)
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/4;
    R = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
    out = mr.rotate3D(R, gy);
    area_y = 0; area_z = 0;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel')
            if out{i}.channel == 'y'
                area_y = out{i}.area;
            elseif out{i}.channel == 'z'
                area_z = out{i}.area;
            end
        end
    end
    testCase.verifyEqual(area_y, gy.area * cos(angle), 'AbsTol', 10);
    testCase.verifyEqual(area_z, gy.area * sin(angle), 'AbsTol', 10);
end

%% Test 45-degree rotation about y splits z into z and x
function test_45deg_y_rotation(testCase)
    gz = mr.makeTrapezoid('z', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/4;
    R = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
    out = mr.rotate3D(R, gz);
    area_z = 0; area_x = 0;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel')
            if out{i}.channel == 'z'
                area_z = out{i}.area;
            elseif out{i}.channel == 'x'
                area_x = out{i}.area;
            end
        end
    end
    testCase.verifyEqual(area_z, gz.area * cos(angle), 'AbsTol', 10);
    testCase.verifyEqual(area_x, gz.area * sin(angle), 'AbsTol', 10);
end

%% Test 70-degree rotation about oblique axis and inverse restore
function test_oblique_rotation_and_inverse(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    % Oblique axis with all components between 0.2 and 0.7
    ax = [0.5, 0.3, 0.6];
    ax = ax / norm(ax); % normalize
    angle = 70 * pi / 180;
    % Rodrigues' rotation formula: R = I*cos(a) + (1-cos(a))*(k*k') + K*sin(a)
    K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
    R = eye(3)*cos(angle) + (1-cos(angle))*(ax'*ax) + K*sin(angle);
    % Forward rotation: should split x-gradient onto all three axes
    out = mr.rotate3D(R, gx);
    area_x = 0; area_y = 0; area_z = 0;
    grads = struct();
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'channel')
            if out{i}.channel == 'x'
                area_x = out{i}.area; grads.x = out{i};
            elseif out{i}.channel == 'y'
                area_y = out{i}.area; grads.y = out{i};
            elseif out{i}.channel == 'z'
                area_z = out{i}.area; grads.z = out{i};
            end
        end
    end
    % Verify all three channels are non-zero (oblique axis splits onto all)
    testCase.verifyGreaterThan(abs(area_x), 1, 'X component should be non-zero');
    testCase.verifyGreaterThan(abs(area_y), 1, 'Y component should be non-zero');
    testCase.verifyGreaterThan(abs(area_z), 1, 'Z component should be non-zero');
    % Verify total area is preserved: sqrt(ax^2+ay^2+az^2) == original
    testCase.verifyEqual(sqrt(area_x^2 + area_y^2 + area_z^2), abs(gx.area), 'AbsTol', 10);
    % Inverse rotation (negative angle about same axis) should restore original
    R_inv = eye(3)*cos(-angle) + (1-cos(-angle))*(ax'*ax) + K*sin(-angle);
    % Build input args: pass all three gradient components
    args = {};
    if isfield(grads, 'x'), args{end+1} = grads.x; end
    if isfield(grads, 'y'), args{end+1} = grads.y; end
    if isfield(grads, 'z'), args{end+1} = grads.z; end
    out2 = mr.rotate3D(R_inv, args{:});
    % After inverse, should be back to x-only
    restored_x = 0; restored_y = 0; restored_z = 0;
    for i = 1:length(out2)
        if isstruct(out2{i}) && isfield(out2{i}, 'channel')
            if out2{i}.channel == 'x'
                restored_x = out2{i}.area;
            elseif out2{i}.channel == 'y'
                restored_y = out2{i}.area;
            elseif out2{i}.channel == 'z'
                restored_z = out2{i}.area;
            end
        end
    end
    testCase.verifyEqual(restored_x, gx.area, 'AbsTol', 15, ...
        'X area should be restored after inverse rotation');
    testCase.verifyEqual(abs(restored_y), 0, 'AbsTol', 15, ...
        'Y area should be ~0 after inverse rotation');
    testCase.verifyEqual(abs(restored_z), 0, 'AbsTol', 15, ...
        'Z area should be ~0 after inverse rotation');
end

%% Test quaternion input
function test_quaternion_input(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/2;
    % Quaternion for 90° about z: [cos(angle/2), 0, 0, sin(angle/2)]
    q = [cos(angle/2), 0, 0, sin(angle/2)];
    out = mr.rotate3D(q, gx);
    testCase.verifyTrue(length(out) >= 1);
end

%% Test non-gradient events pass through
function test_non_grad_passthrough(testCase)
    adc = mr.makeAdc(128, 'Duration', 1e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    R = eye(3);
    out = mr.rotate3D(R, gx, adc);
    % Should have both events
    has_adc = false;
    for i = 1:length(out)
        if isstruct(out{i}) && isfield(out{i}, 'type') && strcmp(out{i}.type, 'adc')
            has_adc = true;
        end
    end
    testCase.verifyTrue(has_adc, 'ADC event should pass through');
end

%% Test duplicate axis gradients throw error
function test_duplicate_axis_error(testCase)
    gx1 = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    gx2 = mr.makeTrapezoid('x', 'Area', 500, 'Duration', 2e-3);
    R = eye(3);
    verifyErrorThrown(testCase, @() mr.rotate3D(R, gx1, gx2));
end

%% Test invalid rotation size throws error
function test_invalid_rotation_error(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    verifyErrorThrown(testCase, @() mr.rotate3D([1 2 3], gx));
end

%% Test with system parameter
function test_with_system(testCase)
    sys = mr.opts();
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    angle = pi/4;
    R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
    out = mr.rotate3D(R, 'system', sys, gx);
    testCase.verifyTrue(~isempty(out));
end

%% Helper
function verifyErrorThrown(testCase, funcHandle)
    didError = false;
    try
        funcHandle();
    catch
        didError = true;
    end
    testCase.verifyTrue(didError, 'Expected an error, but none was thrown.');
end
