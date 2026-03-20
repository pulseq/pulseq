function tests = testSimRf
    tests = functiontests(localfunctions);
end

%% Test 90° block pulse excitation
function test_90deg_excitation(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    [Mz_z, Mz_xy, F] = mr.simRf(rf);
    % On-resonance: find the frequency bin closest to 0
    [~, i0] = min(abs(F));
    % After 90° excitation: Mz ≈ 0 on resonance
    testCase.verifyEqual(Mz_z(i0), 0, 'AbsTol', 0.15, ...
        '90° pulse should tip Mz to ~0 on resonance');
    % Transverse magnetization should be near 1
    testCase.verifyEqual(abs(Mz_xy(i0)), 1, 'AbsTol', 0.15);
end

%% Test 180° block pulse inversion
function test_180deg_inversion(testCase)
    rf = mr.makeBlockPulse(pi, 'Duration', 0.7e-3);
    [Mz_z, ~, F] = mr.simRf(rf);
    [~, i0] = min(abs(F));
    % After 180° inversion: Mz ≈ -1 on resonance
    testCase.verifyEqual(Mz_z(i0), -1, 'AbsTol', 0.15, ...
        '180° pulse should invert Mz to ~-1 on resonance');
end

%% Test output dimensions match
function test_output_dimensions(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    [Mz_z, Mz_xy, F, ref_eff, Mx_xy, My_xy] = mr.simRf(rf);
    n = length(F);
    testCase.verifyEqual(length(Mz_z), n);
    testCase.verifyEqual(length(Mz_xy), n);
    testCase.verifyEqual(length(ref_eff), n);
    testCase.verifyEqual(length(Mx_xy), n);
    testCase.verifyEqual(length(My_xy), n);
end

%% Test sinc pulse profile
function test_sinc_pulse_profile(testCase)
    rf = mr.makeSincPulse(pi/2, 'Duration', 4e-3, 'timeBwProduct', 4);
    [Mz_z, ~, F] = mr.simRf(rf, 0);
    % On-resonance should show excitation
    [~, i0] = min(abs(F));
    testCase.verifyTrue(Mz_z(i0) < 0.3, ...
        'On-resonance longitudinal magnetization should be reduced after 90° excitation');
end
