function tests = testAddRamps
    tests = functiontests(localfunctions);
end

%% Test single-channel trajectory starts and ends at zero
function test_single_channel(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    n = 20;
    % A simple linear k-space trajectory
    k = linspace(0, 100, n);
    kout = mr.addRamps(k, sys);
    % Output should start near 0 and end near 0
    testCase.verifyEqual(kout(1), 0, 'AbsTol', 1e-6);
    testCase.verifyEqual(kout(end), 0, 'AbsTol', 1e-6);
    % Should be longer than input
    testCase.verifyTrue(size(kout, 2) > size(k, 2));
end

%% Test multi-channel cell array
function test_multi_channel(testCase)
    sys = mr.opts();
    n = 20;
    kx = linspace(0, 100, n);
    ky = linspace(0, 50, n);
    [kx_out, ky_out] = mr.addRamps({kx, ky}, sys);
    % Both should start and end at 0
    testCase.verifyEqual(kx_out(1), 0, 'AbsTol', 1e-3);
    testCase.verifyEqual(kx_out(end), 0, 'AbsTol', 1e-3);
    testCase.verifyEqual(ky_out(1), 0, 'AbsTol', 1e-3);
    testCase.verifyEqual(ky_out(end), 0, 'AbsTol', 1e-3);
    % Same length
    testCase.verifyEqual(length(kx_out), length(ky_out));
end

%% Test with RF padding
function test_rf_padding(testCase)
    sys = mr.opts();
    n = 20;
    k = linspace(0, 100, n);
    rf_shape = ones(1, n * 10); % RF sampled at higher rate
    [kout, rf_out] = mr.addRamps(k, sys, 'rf', rf_shape);
    % RF output should be longer (zero-padded for ramps)
    testCase.verifyTrue(length(rf_out) > length(rf_shape));
    % First and last samples of RF should be zero (ramp padding)
    testCase.verifyEqual(rf_out(1), 0, 'AbsTol', 1e-10);
    testCase.verifyEqual(rf_out(end), 0, 'AbsTol', 1e-10);
end

%% Test trajectory that already starts/ends at zero
function test_zero_start_end(testCase)
    sys = mr.opts();
    n = 20;
    k = [zeros(1,3) linspace(0, 100, n) zeros(1,3)];
    kout = mr.addRamps(k, sys);
    testCase.verifyEqual(kout(1), 0, 'AbsTol', 1e-6);
    testCase.verifyEqual(kout(end), 0, 'AbsTol', 1e-6);
end
