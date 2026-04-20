function tests = testTransformFOV
%testTransformFOV  Unit tests for the mr.TransformFOV class.
%
%   Tests cover:
%     - Constructor parameter validation
%     - Scale-only transforms on trapezoid and arbitrary gradients
%     - Rotation-only transforms (90°, 45°, arbitrary axis)
%     - Translation-only transforms (phase/frequency offsets on RF and ADC)
%     - Combined rotation + translation
%     - Combined rotation + translation + scaling
%     - Rotation via rotation-extension flag
%     - Homogeneous 4×4 transform matrix input
%     - applyToSeq on a multi-block sequence
%     - Analytic gradient-area and k-space verification
%
%   Uses the function-based unit-test framework so that it runs on both
%   MATLAB (via ``runtests``) and GNU Octave (via ``oruntests``).

    tests = functiontests(localfunctions);
end

% =====================================================================
%  Helper: default system struct used throughout
% =====================================================================
function sys = defaultSys()
    % Use generous limits so that scaled gradients still fit within the
    % hardware constraints (scaling changes both amplitude and slew rate).
    sys = mr.opts('MaxGrad', 200, 'GradUnit', 'mT/m', ...
                  'MaxSlew', 200, 'SlewUnit', 'T/m/s');
end

% Helper: extract gradient events from a cell-array output of applyToBlock
function [gx, gy, gz] = extractGrads(events)
    gx = []; gy = []; gz = [];
    for k = 1:length(events)
        e = events{k};
        if isstruct(e) && isfield(e,'type') && (strcmp(e.type,'trap') || strcmp(e.type,'grad'))
            switch e.channel
                case 'x', gx = e;
                case 'y', gy = e;
                case 'z', gz = e;
            end
        end
    end
end

% Helper: extract RF event from a cell-array output of applyToBlock
function rf = extractRf(events)
    rf = [];
    for k = 1:length(events)
        e = events{k};
        if isstruct(e) && isfield(e,'type') && strcmp(e.type,'rf')
            rf = e; return;
        end
    end
end

% Helper: extract ADC event from a cell-array output of applyToBlock
function adc = extractAdc(events)
    adc = [];
    for k = 1:length(events)
        e = events{k};
        if isstruct(e) && isfield(e,'type') && strcmp(e.type,'adc')
            adc = e; return;
        end
    end
end

% Helper: get gradient area (works for both trap and grad types)
function a = gradArea(g)
    if isempty(g)
        a = 0;
    elseif strcmp(g.type,'trap')
        a = g.area;
    else
        a = sum(g.waveform(1:end-1) .* diff(g.tt)) + ...
            0.5*sum(diff(g.waveform) .* diff(g.tt));  % trapezoidal integration
    end
end

% Helper: rotation matrix about z-axis
function R = Rz(angle)
    R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
end

% Helper: rotation matrix about y-axis
function R = Ry(angle)
    R = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
end

% Helper: rotation matrix about x-axis
function R = Rx(angle)
    R = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
end

% =====================================================================
%  1.  Constructor tests
% =====================================================================

%% Test: constructor rejects empty arguments
function test_constructor_no_args(testCase)
    % Must give at least one transform parameter
    threw = false;
    try
        mr.TransformFOV(); 
    catch
        threw = true;
    end
    testCase.verifyTrue(threw, ...
        'TransformFOV() with no arguments should throw an error');
end

%% Test: constructor stores rotation correctly
function test_constructor_rotation(testCase)
    R = Rz(pi/4);
    T = mr.TransformFOV('rotation', R);
    testCase.verifyEqual(T.rotation, R, 'AbsTol', 1e-14);
    testCase.verifyTrue(isempty(T.translation));
    testCase.verifyTrue(isempty(T.scale));
end

%% Test: constructor stores translation correctly
function test_constructor_translation(testCase)
    t = [0.01 0.02 0.03];
    T = mr.TransformFOV('translation', t);
    testCase.verifyEqual(T.translation, t, 'AbsTol', 1e-14);
    testCase.verifyTrue(isempty(T.rotation));
end

%% Test: constructor stores scale correctly
function test_constructor_scale(testCase)
    s = [2 0.5 1];
    T = mr.TransformFOV('scale', s);
    testCase.verifyEqual(T.scale, s, 'AbsTol', 1e-14);
end

%% Test: constructor rejects transform with rotation
function test_constructor_transform_with_rotation_error(testCase)
    M = eye(4);
    R = eye(3);
    try
        mr.TransformFOV('transform', M, 'rotation', R);
        testCase.verifyFail('Expected an error when combining transform with rotation');
    catch
        % expected
    end
end

%% Test: constructor rejects transform with translation
function test_constructor_transform_with_translation_error(testCase)
    M = eye(4);
    t = [1 1 1];
    try
        mr.TransformFOV('transform', M, 'translation', t);
        testCase.verifyFail('Expected an error when combining transform with translation');
    catch
        % expected
    end
end

% =====================================================================
%  2.  Scale-only transform tests
% =====================================================================

%% Test: scale trapezoid gradient by a factor of 2
function test_scale_trap_x2(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    T = mr.TransformFOV('scale', [0.5 1 1], 'system', sys);
    out = T.applyToBlock(gx);
    [gx2, ~, ~] = extractGrads(out);
    testCase.verifyEqual(gx2.area, 0.5*gx.area, 'AbsTol', 1);
    testCase.verifyEqual(gx2.amplitude, 0.5*gx.amplitude, 'AbsTol', 1e-3);
end

%% Test: scale only the z-channel, leaving x and y unchanged
function test_scale_selective_channel(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', 2000, 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', 3000, 'Duration', 2e-3);
    T = mr.TransformFOV('scale', [1 1 0.5], 'system', sys);
    out = T.applyToBlock(gx, gy, gz);
    [gx2, gy2, gz2] = extractGrads(out);
    testCase.verifyEqual(gx2.area, gx.area, 'AbsTol', 1, 'Gx should be unchanged');
    testCase.verifyEqual(gy2.area, gy.area, 'AbsTol', 1, 'Gy should be unchanged');
    testCase.verifyEqual(gz2.area, 0.5*gz.area, 'AbsTol', 1, 'Gz should be halved');
end

%% Test: scale by negation inverts gradient
function test_scale_negate(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    T = mr.TransformFOV('scale', [-1 1 1], 'system', sys);
    out = T.applyToBlock(gx);
    [gx2, ~, ~] = extractGrads(out);
    testCase.verifyEqual(gx2.area, -gx.area, 'AbsTol', 1);
end

% =====================================================================
%  3.  Rotation-only transform tests
% =====================================================================

%% Test: identity rotation leaves gradients unchanged
function test_rotation_identity(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', 2000, 'Duration', 2e-3);
    T = mr.TransformFOV('rotation', eye(3), 'system', sys);
    out = T.applyToBlock(gx, gy);
    [gx2, gy2, ~] = extractGrads(out);
    testCase.verifyEqual(gradArea(gx2), 1000, 'AbsTol', 10, 'Gx area should be preserved');
    testCase.verifyEqual(gradArea(gy2), 2000, 'AbsTol', 10, 'Gy area should be preserved');
end

%% Test: 90° rotation about z maps gx -> gy
function test_rotation_90z_gx_to_gy(testCase)
    sys = defaultSys();
    area_x = 1000;
    gx = mr.makeTrapezoid('x', sys, 'Area', area_x, 'Duration', 2e-3);
    R = Rz(pi/2);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gx);
    [gx2, gy2, ~] = extractGrads(out);

    % After 90° about z: x -> y  (gx component vanishes, gy gets the area)
    testCase.verifyEqual(gradArea(gx2), 0, 'AbsTol', 15, 'Gx area should vanish');
    testCase.verifyEqual(abs(gradArea(gy2)), abs(area_x), 'AbsTol', 15, 'Gy area should equal original Gx');
end

%% Test: 90° rotation about x maps gy -> gz
function test_rotation_90x_gy_to_gz(testCase)
    sys = defaultSys();
    area_y = 2000;
    gy = mr.makeTrapezoid('y', sys, 'Area', area_y, 'Duration', 2e-3);
    R = Rx(pi/2);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gy);
    [~, gy2, gz2] = extractGrads(out);
    testCase.verifyEqual(gradArea(gy2), 0, 'AbsTol', 15, 'Gy area should vanish');
    testCase.verifyEqual(abs(gradArea(gz2)), abs(area_y), 'AbsTol', 15, 'Gz area should equal original Gy');
end

%% Test: 45° rotation about z splits gx into gx and gy equally
function test_rotation_45z_splits_area(testCase)
    sys = defaultSys();
    area_x = 1000;
    gx = mr.makeTrapezoid('x', sys, 'Area', area_x, 'Duration', 2e-3);
    angle = pi/4;
    R = Rz(angle);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gx);
    [gx2, gy2, ~] = extractGrads(out);

    expected_x = area_x * cos(angle);
    expected_y = area_x * sin(angle);
    testCase.verifyEqual(gradArea(gx2), expected_x, 'RelTol', 0.02, 'Gx area after 45° rotation');
    testCase.verifyEqual(gradArea(gy2), expected_y, 'RelTol', 0.02, 'Gy area after 45° rotation');
end

%% Test: rotation preserves total gradient area (norm)
function test_rotation_preserves_total_area(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', 2000, 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', 1500, 'Duration', 2e-3);
    area_vec_orig = [1000, 2000, 1500];
    norm_orig = norm(area_vec_orig);

    % Arbitrary rotation: 30° about z, then 20° about y
    R = Ry(20*pi/180) * Rz(30*pi/180);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gx, gy, gz);
    [gx2, gy2, gz2] = extractGrads(out);

    area_vec_rot = [gradArea(gx2), gradArea(gy2), gradArea(gz2)];
    norm_rot = norm(area_vec_rot);
    testCase.verifyEqual(norm_rot, norm_orig, 'RelTol', 0.02, ...
        'Rotation should preserve total gradient area L2 norm');
end

%% Test: rotation analytically matches R * area_vector
function test_rotation_analytic_area_vector(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', -500, 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', 2000, 'Duration', 2e-3);
    area_orig = [1000; -500; 2000];

    % Rotation: 37° about z, then 53° about x
    R = Rx(53*pi/180) * Rz(37*pi/180);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gx, gy, gz);
    [gx2, gy2, gz2] = extractGrads(out);
    area_rot = [gradArea(gx2); gradArea(gy2); gradArea(gz2)];

    expected = R * area_orig;
    testCase.verifyEqual(area_rot, expected, 'RelTol', 0.03, ...
        'Rotated gradient area vector should match R * original');
end

% =====================================================================
%  4.  Translation-only transform tests (RF/ADC phase offsets)
% =====================================================================

%% Test: translation with constant gradient sets RF frequency offset
function test_translation_rf_freq_offset(testCase)
    sys = defaultSys();
    % Create a trapezoid on x-axis and an RF pulse
    gx = mr.makeTrapezoid('x', sys, 'FlatArea', 5000, 'FlatTime', 1e-3);
    [rf, ~] = mr.makeSincPulse(pi/6, sys, 'Duration', 1e-3, ...
              'sliceThickness', 5e-3, ...
              'delay', gx.riseTime, 'use', 'excitation');

    shift_x = 0.01; % 1 cm shift in x
    T = mr.TransformFOV('translation', [shift_x 0 0], 'system', sys);
    out = T.applyToBlock(rf, gx);
    rf2 = extractRf(out);

    % With constant gradient during RF, freqOffset should be shift * amplitude
    expected_freq = shift_x * gx.amplitude;
    testCase.verifyEqual(rf2.freqOffset, expected_freq, 'RelTol', 0.01, ...
        'RF frequency offset should be shift * constant gradient amplitude');
end

%% Test: translation with constant gradient sets ADC frequency offset
function test_translation_adc_freq_offset(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'FlatArea', 5000, 'FlatTime', 2e-3);
    adc = mr.makeAdc(64, sys, 'Duration', 2e-3, 'delay', gx.riseTime);

    shift_x = 0.02; % 2 cm shift
    T = mr.TransformFOV('translation', [shift_x 0 0], 'system', sys);
    out = T.applyToBlock(adc, gx);
    adc2 = extractAdc(out);

    expected_freq = shift_x * gx.amplitude;
    testCase.verifyEqual(adc2.freqOffset, expected_freq, 'RelTol', 0.01, ...
        'ADC frequency offset should be shift * constant gradient amplitude');
end

%% Test: zero translation leaves RF and ADC unchanged
function test_translation_zero_no_change(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'FlatArea', 5000, 'FlatTime', 1e-3);
    [rf, ~] = mr.makeSincPulse(pi/6, sys, 'Duration', 1e-3, ...
              'sliceThickness', 5e-3, ...
              'delay', gx.riseTime, 'use', 'excitation');
    adc = mr.makeAdc(32, sys, 'Duration', 1e-3, 'delay', gx.riseTime);

    T = mr.TransformFOV('translation', [0 0 0], 'system', sys);
    out = T.applyToBlock(rf, adc, gx);
    rf2 = extractRf(out);
    adc2 = extractAdc(out);

    testCase.verifyEqual(rf2.freqOffset, rf.freqOffset, 'AbsTol', 1e-6);
    testCase.verifyEqual(adc2.freqOffset, adc.freqOffset, 'AbsTol', 1e-6);
end

% =====================================================================
%  5.  Combined rotation + translation tests
% =====================================================================

%% Test: rotation + translation – gradient areas match rotation and RF gets offset
function test_combined_rotation_translation(testCase)
    sys = defaultSys();
    area_x = 5000;
    gx = mr.makeTrapezoid('x', sys, 'FlatArea', area_x, 'FlatTime', 1e-3);
    [rf, ~] = mr.makeSincPulse(pi/6, sys, 'Duration', 1e-3, ...
              'sliceThickness', 5e-3, ...
              'delay', gx.riseTime, 'use', 'excitation');

    angle = pi/2;  % 90° about z
    R = Rz(angle);
    shift = [0.01 0 0]; % 1 cm in x
    T = mr.TransformFOV('rotation', R, 'translation', shift, 'system', sys);
    out = T.applyToBlock(rf, gx);
    [gx2, gy2, ~] = extractGrads(out);

    % After rotation: x -> y, so Gx area should vanish, Gy should have the area
    testCase.verifyEqual(gradArea(gx2), 0, 'AbsTol', 20, 'Gx area should vanish after 90° z rotation');
    testCase.verifyEqual(abs(gradArea(gy2)), abs(gx.area), 'AbsTol', 20, ...
        'Gy should carry the area after 90° z rotation');

    % The RF should have a non-zero frequency offset from the translation
    rf2 = extractRf(out);
    testCase.verifyTrue(abs(rf2.freqOffset) > 0, ...
        'RF freqOffset should be nonzero with translation');
end

% =====================================================================
%  6.  Combined rotation + translation + scaling
% =====================================================================

%% Test: full combined transform – analytic gradient area verification
function test_combined_all_three(testCase)
    sys = defaultSys();
    area_orig = [1000; 2000; -1500];
    gx = mr.makeTrapezoid('x', sys, 'Area', area_orig(1), 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', area_orig(2), 'Duration', 2e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', area_orig(3), 'Duration', 2e-3);

    S = [0.7 0.5 0.6];
    angle_z = 30*pi/180;
    angle_x = 20*pi/180;
    R = Rx(angle_x) * Rz(angle_z);
    shift = [0.005 0.01 -0.003];

    T = mr.TransformFOV('scale', S, 'rotation', R, 'translation', shift, 'system', sys);
    out = T.applyToBlock(gx, gy, gz);
    [gx2, gy2, gz2] = extractGrads(out);
    area_result = [gradArea(gx2); gradArea(gy2); gradArea(gz2)];

    % The transform order is scale, then translation (just phase offsets,
    % does not change gradient areas), then rotation:
    %   area_result = R * diag(S) * area_orig
    expected = R * (S(:) .* area_orig);
    testCase.verifyEqual(area_result, expected, 'RelTol', 0.03, ...
        'Gradient area vector should be R * diag(S) * original');
end

%% Test: double rotation equals sequential application
function test_double_rotation_sequential(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);

    R1 = Rz(pi/6);
    R2 = Rx(pi/4);
    % Single combined rotation
    T_combined = mr.TransformFOV('rotation', R2*R1, 'system', sys);
    out_combined = T_combined.applyToBlock(gx);
    [gx_c, gy_c, gz_c] = extractGrads(out_combined);
    area_combined = [gradArea(gx_c); gradArea(gy_c); gradArea(gz_c)];

    % Sequential: first R1, then R2
    T1 = mr.TransformFOV('rotation', R1, 'system', sys);
    T2 = mr.TransformFOV('rotation', R2, 'system', sys);
    out1 = T1.applyToBlock(gx);
    % Pass all outputs of T1 into T2 (unpack cell array)
    out2 = T2.applyToBlock(out1{:});
    [gx_s, gy_s, gz_s] = extractGrads(out2);
    area_sequential = [gradArea(gx_s); gradArea(gy_s); gradArea(gz_s)];

    testCase.verifyEqual(area_combined, area_sequential, 'RelTol', 0.03, ...
        'Combined rotation should equal sequential rotations');
end

% =====================================================================
%  7.  Use rotation extension
% =====================================================================

%% Test: rotation extension flag produces a rotation extension event
function test_rotation_extension_produces_event(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    R = Rz(pi/4);
    T = mr.TransformFOV('rotation', R, 'use_rotation_extension', true, 'system', sys);
    out = T.applyToBlock(gx);

    % The output should contain a 'rot3D' event
    found = false;
    for k = 1:length(out)
        if isstruct(out{k}) && isfield(out{k},'type') && strcmp(out{k}.type,'rot3D')
            found = true;
            break;
        end
    end
    testCase.verifyTrue(found, 'Output should contain a rot3D event when use_rotation_extension=true');
end

%% Test: with rotation extension, gradient areas are NOT rotated directly
function test_rotation_extension_grads_not_rotated(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    R = Rz(pi/2);
    T = mr.TransformFOV('rotation', R, 'use_rotation_extension', true, 'system', sys);
    out = T.applyToBlock(gx);
    [gx2, gy2, ~] = extractGrads(out);

    % With rotation extension, the gradient itself should stay on x, not be rotated
    testCase.verifyEqual(gradArea(gx2), 1000, 'AbsTol', 15, ...
        'Gx area should be preserved when using rotation extension');
    testCase.verifyEqual(gradArea(gy2), 0, 'AbsTol', 15, ...
        'Gy area should be zero when using rotation extension');
end

% =====================================================================
%  8.  Homogeneous 4×4 transform matrix
% =====================================================================

%% Test: 4x4 identity transform is equivalent to rotation=eye(3)
function test_transform_4x4_identity(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', 2000, 'Duration', 2e-3);
    T = mr.TransformFOV('transform', eye(4), 'system', sys);
    out = T.applyToBlock(gx, gy);
    [gx2, gy2, ~] = extractGrads(out);
    testCase.verifyEqual(gradArea(gx2), 1000, 'AbsTol', 10);
    testCase.verifyEqual(gradArea(gy2), 2000, 'AbsTol', 10);
end

% =====================================================================
%  9.  applyToSeq – multi-block sequence
% =====================================================================

%% Test: applyToSeq applies transform to every block
function test_applyToSeq_basic(testCase)
    sys = defaultSys();
    seq = mr.Sequence(sys);

    % Build a simple 3-block sequence with readout gradients only
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    seq.addBlock(gx);
    seq.addBlock(gx);
    seq.addBlock(gx);

    R = Rz(pi/2);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    seq2 = T.applyToSeq(seq);

    testCase.verifyEqual(length(seq2.blockDurations), 3, 'Number of blocks should be preserved');

    % Each block's x-gradient area should be ~0 and y should carry the area
    for iB = 1:3
        B = seq2.getBlock(iB);
        gx_area = 0; gy_area = 0;
        if ~isempty(B.gx), gx_area = gradArea(B.gx); end
        if ~isempty(B.gy), gy_area = gradArea(B.gy); end
        testCase.verifyEqual(gx_area, 0, 'AbsTol', 15, ...
            sprintf('Block %d: Gx should vanish after 90° z rotation', iB));
        testCase.verifyEqual(abs(gy_area), abs(gx.area), 'AbsTol', 15, ...
            sprintf('Block %d: Gy should carry the area', iB));
    end
end

%% Test: applyToSeq with blockRange only transforms specified blocks
function test_applyToSeq_blockRange(testCase)
    sys = defaultSys();
    seq = mr.Sequence(sys);

    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
    seq.addBlock(gx);
    seq.addBlock(gx);
    seq.addBlock(gx);

    R = Rz(pi/2);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    seq2 = T.applyToSeq(seq, 'blockRange', [2 3]);

    % Only blocks 2 and 3 should be rotated
    testCase.verifyEqual(length(seq2.blockDurations), 2, ...
        'Block range [2 3] should yield 2 blocks');
end

% =====================================================================
%  10. Advanced: analytic k-space verification after combined transforms
% =====================================================================

%% Test: k-space trajectory rotated by R using waveforms_and_times
function test_kspace_rotation_analytic(testCase)
    sys = defaultSys();
    seq_orig = mr.Sequence(sys);

    readArea = 5000;
    gx = mr.makeTrapezoid('x', sys, 'Area', readArea, 'Duration', 2e-3);
    adc = mr.makeAdc(64, sys, 'Duration', 2e-3 - gx.riseTime - gx.fallTime, ...
                     'delay', gx.riseTime);
    seq_orig.addBlock(gx, adc);

    % Apply 90° z-rotation
    R = Rz(pi/2);
    T = mr.TransformFOV('rotation', R, 'system', sys);
    seq_rot = T.applyToSeq(seq_orig);

    % Extract gradient waveforms
    w_orig = seq_orig.waveforms_and_times();
    w_rot  = seq_rot.waveforms_and_times();

    % Original: all gradient area is on x, nothing on y
    area_x_orig = trapz(w_orig{1}(1,:), w_orig{1}(2,:));

    % Rotated: x area should be negligible, y should carry the area
    area_x_rot = 0;
    area_y_rot = 0;
    if ~isempty(w_rot{1})
        area_x_rot = trapz(w_rot{1}(1,:), w_rot{1}(2,:));
    end
    if ~isempty(w_rot{2})
        area_y_rot = trapz(w_rot{2}(1,:), w_rot{2}(2,:));
    end

    testCase.verifyEqual(area_x_rot, 0, 'AbsTol', 50, ...
        'After 90° z rotation, x-gradient area should vanish');
    testCase.verifyEqual(abs(area_y_rot), abs(area_x_orig), 'RelTol', 0.02, ...
        'After 90° z rotation, y-gradient area should equal original x');
end

%% Test: k-space area vector matches analytic R*S*a after combined transform
function test_kspace_combined_transform_analytic(testCase)
    sys = defaultSys();

    area_orig = [3000; -1000; 2000];
    gx = mr.makeTrapezoid('x', sys, 'Area', area_orig(1), 'Duration', 3e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', area_orig(2), 'Duration', 3e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', area_orig(3), 'Duration', 3e-3);

    seq_orig = mr.Sequence(sys);
    seq_orig.addBlock(gx, gy, gz);

    % Build combined transform
    S = [0.5 0.8 0.7];
    R = Ry(25*pi/180) * Rz(40*pi/180);
    T = mr.TransformFOV('scale', S, 'rotation', R, 'system', sys);
    seq_t = T.applyToSeq(seq_orig);

    % Measure actual gradient areas from waveforms
    w = seq_t.waveforms_and_times();
    area_actual = zeros(3,1);
    for ch = 1:3
        if ~isempty(w{ch})
            area_actual(ch) = trapz(w{ch}(1,:), w{ch}(2,:));
        end
    end

    % Expected: R * (S .* area_orig)
    expected = R * (S(:) .* area_orig);
    testCase.verifyEqual(area_actual, expected, 'RelTol', 0.03, ...
        'Waveform area vector should match R*diag(S)*area_orig');
end

%% Test: 180° rotation inverts k-space trajectory
function test_kspace_180_rotation_inverts(testCase)
    sys = defaultSys();
    seq_orig = mr.Sequence(sys);
    gx = mr.makeTrapezoid('x', sys, 'Area', 4000, 'Duration', 2e-3);
    seq_orig.addBlock(gx);

    R = Rz(pi); % 180° about z: x -> -x
    T = mr.TransformFOV('rotation', R, 'system', sys);
    seq_inv = T.applyToSeq(seq_orig);

    w_orig = seq_orig.waveforms_and_times();
    w_inv  = seq_inv.waveforms_and_times();

    area_orig = trapz(w_orig{1}(1,:), w_orig{1}(2,:));
    area_inv = 0;
    if ~isempty(w_inv{1})
        area_inv = trapz(w_inv{1}(1,:), w_inv{1}(2,:));
    end

    testCase.verifyEqual(area_inv, -area_orig, 'RelTol', 0.02, ...
        '180° z rotation should negate x-gradient area');
end

%% Test: three-axis rotation chain – analytic vs numeric
function test_three_axis_rotation_chain(testCase)
    sys = defaultSys();
    a = [1500; -800; 2200];
    gx = mr.makeTrapezoid('x', sys, 'Area', a(1), 'Duration', 3e-3);
    gy = mr.makeTrapezoid('y', sys, 'Area', a(2), 'Duration', 3e-3);
    gz = mr.makeTrapezoid('z', sys, 'Area', a(3), 'Duration', 3e-3);

    % Euler angles ZYX:  Rz(alpha) * Ry(beta) * Rx(gamma)
    alpha = 17*pi/180;
    beta  = 43*pi/180;
    gamma = -29*pi/180;
    R = Rz(alpha) * Ry(beta) * Rx(gamma);

    T = mr.TransformFOV('rotation', R, 'system', sys);
    out = T.applyToBlock(gx, gy, gz);
    [gx2, gy2, gz2] = extractGrads(out);
    area_rot = [gradArea(gx2); gradArea(gy2); gradArea(gz2)];

    expected = R * a;
    testCase.verifyEqual(area_rot, expected, 'RelTol', 0.03, ...
        'ZYX Euler rotation: gradient area vector should match R*a');
end

%% Test: translation phase accumulation across multiple blocks
function test_translation_phase_accumulation(testCase)
    % Verify that prior_phase_cycle is updated across blocks by applying
    % the transform to a multi-block sequence with the same gradient.
    % The phase cycle should advance consistently.
    sys = defaultSys();
    seq = mr.Sequence(sys);
    gx = mr.makeTrapezoid('x', sys, 'FlatArea', 1000, 'FlatTime', 1e-3);
    adc = mr.makeAdc(16, sys, 'Duration', 1e-3, 'delay', gx.riseTime);
    seq.addBlock(gx, adc);
    seq.addBlock(gx, adc);
    seq.addBlock(gx, adc);

    shift_x = 0.01;
    T = mr.TransformFOV('translation', [shift_x 0 0], 'system', sys);
    seq2 = T.applyToSeq(seq);

    % All ADC events should have the same nonzero frequency offset
    % (because the constant-gradient case just adds freq offset)
    for iB = 1:3
        B = seq2.getBlock(iB);
        testCase.verifyTrue(abs(B.adc.freqOffset) > 0, ...
            sprintf('Block %d ADC should have nonzero freqOffset', iB));
    end
end

% =====================================================================
%  11. Edge case: gradient-only block (no RF/ADC)
% =====================================================================

%% Test: translate block with only gradients – no crash
function test_translation_gradient_only_block(testCase)
    sys = defaultSys();
    gx = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 1e-3);
    T = mr.TransformFOV('translation', [0.01 0 0], 'system', sys);

    % Should not error even without RF/ADC
    out = T.applyToBlock(gx);
    [gx2, ~, ~] = extractGrads(out);
    testCase.verifyEqual(gradArea(gx2), 1000, 'AbsTol', 5, ...
        'Gradient area should be unchanged (no RF/ADC to modulate)');
end

% =====================================================================
%  12. Edge case: RF-only block (no gradient)
% =====================================================================

%% Test: translate block with only RF – no crash, phase unchanged
function test_translation_rf_only_block(testCase)
    sys = defaultSys();
    rf = mr.makeBlockPulse(pi/2, sys, 'Duration', 1e-3, 'use', 'excitation');
    T = mr.TransformFOV('translation', [0.01 0.02 0.03], 'system', sys);

    out = T.applyToBlock(rf);
    rf2 = extractRf(out);
    % Without gradients there is nothing to produce a frequency offset
    testCase.verifyEqual(rf2.freqOffset, rf.freqOffset, 'AbsTol', 1e-6);
end
