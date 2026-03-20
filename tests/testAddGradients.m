function tests = testAddGradients
    tests = functiontests(localfunctions);
end

%% Test sum of two identical trapezoids doubles amplitude
function test_sum_identical_traps(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    gsum = mr.addGradients({g, g});
    testCase.verifyEqual(gsum.amplitude, 2*g.amplitude, 'AbsTol', 1e-6);
    testCase.verifyEqual(gsum.area, 2*g.area, 'AbsTol', 1e-3);
    testCase.verifyEqual(gsum.flatArea, 2*g.flatArea, 'AbsTol', 1e-3);
end

%% Test opposing trapezoids cancel
function test_opposing_traps(testCase)
    g1 = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g2 = mr.makeTrapezoid('x', 'Area', -1000, 'Duration', 1e-3);
    gsum = mr.addGradients({g1, g2});
    testCase.verifyEqual(gsum.area, 0, 'AbsTol', 1);
end

%% Test non-cell input error
function test_non_cell_error(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    testCase.verifyError(@() mr.addGradients(g), ?MException);
end

%% Test single gradient error
function test_single_gradient_error(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    testCase.verifyError(@() mr.addGradients({g}), ?MException);
end

%% Test different channels error
function test_different_channels_error(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 1e-3);
    testCase.verifyError(@() mr.addGradients({gx, gy}), ?MException);
end

%% Test trap + extended trap
function test_trap_plus_extended(testCase)
    sys = mr.opts();
    g_trap = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g_ext = mr.makeExtendedTrapezoid('x', ...
        'Times', [0 2e-4 8e-4 1e-3], ...
        'Amplitudes', [0 50000 50000 0]);
    gsum = mr.addGradients({g_trap, g_ext}, sys);
    % Result should be an extended trapezoid
    testCase.verifyTrue(strcmp(gsum.type, 'grad'));
    % Duration should be the max of both
    dur_sum = mr.calcDuration(gsum);
    dur_max = max(mr.calcDuration(g_trap), mr.calcDuration(g_ext));
    testCase.verifyEqual(dur_sum, dur_max, 'AbsTol', sys.gradRasterTime);
end

%% Test three trapezoids
function test_three_traps(testCase)
    g1 = mr.makeTrapezoid('y', 'Area', 500, 'Duration', 1e-3);
    g2 = mr.makeTrapezoid('y', 'Area', 300, 'Duration', 1e-3);
    g3 = mr.makeTrapezoid('y', 'Area', 200, 'Duration', 1e-3);
    gsum = mr.addGradients({g1, g2, g3});
    testCase.verifyEqual(gsum.area, g1.area + g2.area + g3.area, 'AbsTol', 1);
end
