function tests = testTraj2grad
    tests = functiontests(localfunctions);
end

%% Test linear trajectory gives constant gradient
function test_linear_trajectory(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    n = 20;
    k = linspace(0, 1000, n);  % linear k-space
    [g, sr] = mr.traj2grad(k, 'system', sys);
    % Gradient should be constant: dk/dt = 1000/(n-1) / grt
    expected_g = (k(2) - k(1)) / grt;
    testCase.verifyEqual(g(1, :), expected_g * ones(1, size(g, 2)), 'RelTol', 0.01);
end

%% Test output dimensions
function test_output_dimensions(testCase)
    sys = mr.opts();
    n = 30;
    k = linspace(0, 500, n);
    [g, sr] = mr.traj2grad(k, 'system', sys);
    testCase.verifyEqual(size(g, 2), n-1, 'Gradient should have n-1 samples');
    testCase.verifyEqual(size(sr, 2), n-1, 'Slew rate should have n-1 samples');
end

%% Test multi-channel trajectory
function test_multi_channel(testCase)
    sys = mr.opts();
    n = 20;
    kx = linspace(0, 500, n);
    ky = linspace(0, 300, n);
    k = [kx; ky];
    [g, sr] = mr.traj2grad(k, 'system', sys);
    testCase.verifyEqual(size(g, 1), 2, 'Should have 2 gradient channels');
    testCase.verifyEqual(size(sr, 1), 2, 'Should have 2 slew rate channels');
end

%% Test zero trajectory gives zero gradient
function test_zero_trajectory(testCase)
    sys = mr.opts();
    k = zeros(1, 10);
    [g, ~] = mr.traj2grad(k, 'system', sys);
    testCase.verifyEqual(g, zeros(size(g)), 'AbsTol', 1e-10);
end

%% Test custom RasterTime
function test_custom_raster(testCase)
    sys = mr.opts();
    custom_raster = 5e-6;
    k = linspace(0, 100, 20);
    [g1, ~] = mr.traj2grad(k, 'system', sys, 'RasterTime', custom_raster);
    [g2, ~] = mr.traj2grad(k, 'system', sys);
    % Different raster times should give different gradient magnitudes
    % g = dk / raster, so g1 = 2*g2 (since custom is half of default)
    testCase.verifyEqual(g1(1), g2(1) * sys.gradRasterTime / custom_raster, 'RelTol', 0.01);
end
