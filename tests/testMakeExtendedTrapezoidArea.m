function tests = testMakeExtendedTrapezoidArea
    tests = functiontests(localfunctions);
end

%% Test small area from zero
function test_small_area(testCase)
    sys = mr.opts();
    A = 100;  % small gradient area
    [grad, times, amplitudes] = mr.makeExtendedTrapezoidArea('x', 0, 0, A, sys);
    testCase.verifyEqual(grad.type, 'grad');
    testCase.verifyEqual(grad.channel, 'x');
    testCase.verifyEqual(grad.area, A, 'AbsTol', 1);
end

%% Test large area (needs flat top)
function test_large_area(testCase)
    sys = mr.opts();
    A = 50000;  % large area requiring flat top
    [grad, times, amplitudes] = mr.makeExtendedTrapezoidArea('x', 0, 0, A, sys);
    testCase.verifyEqual(grad.area, A, 'AbsTol', 1);
    % Should have 4 time points (ramp-up, flat, ramp-down)
    testCase.verifyGreaterThanOrEqual(length(times), 3);
end

%% Test with non-zero start/end gradients
function test_nonzero_gs_ge(testCase)
    sys = mr.opts();
    Gs = 10000;  % starting gradient
    Ge = 5000;   % ending gradient
    A = 2000;
    [grad, ~, ~] = mr.makeExtendedTrapezoidArea('x', Gs, Ge, A, sys);
    testCase.verifyEqual(grad.area, A, 'AbsTol', 1);
end

%% Test negative area
function test_negative_area(testCase)
    sys = mr.opts();
    A = -5000;
    [grad, ~, ~] = mr.makeExtendedTrapezoidArea('x', 0, 0, A, sys);
    testCase.verifyEqual(grad.area, A, 'AbsTol', 1);
end

%% Test gradient amplitude within limits
function test_amplitude_within_limits(testCase)
    sys = mr.opts();
    A = 10000;
    [grad, ~, amplitudes] = mr.makeExtendedTrapezoidArea('x', 0, 0, A, sys);
    testCase.verifyTrue(all(abs(amplitudes) <= sys.maxGrad), ...
        'Gradient amplitudes should be within system limits');
end

%% Test different channels
function test_channel_y(testCase)
    sys = mr.opts();
    [grad, ~, ~] = mr.makeExtendedTrapezoidArea('y', 0, 0, 500, sys);
    testCase.verifyEqual(grad.channel, 'y');
end
