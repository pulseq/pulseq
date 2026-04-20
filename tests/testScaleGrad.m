function tests = testScaleGrad
    tests = functiontests(localfunctions);
end

%% Test scale trapezoid by 2
function test_scale_trap_by_2(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g2 = mr.scaleGrad(g, 2);
    testCase.verifyEqual(g2.amplitude, 2*g.amplitude, 'AbsTol', 1e-6);
    testCase.verifyEqual(g2.area, 2*g.area, 'AbsTol', 1e-3);
    testCase.verifyEqual(g2.flatArea, 2*g.flatArea, 'AbsTol', 1e-3);
end

%% Test scale trapezoid by -1 (negate)
function test_scale_trap_negate(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g_neg = mr.scaleGrad(g, -1);
    testCase.verifyEqual(g_neg.amplitude, -g.amplitude, 'AbsTol', 1e-6);
    testCase.verifyEqual(g_neg.area, -g.area, 'AbsTol', 1e-3);
end

%% Test scale by 0
function test_scale_by_zero(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g0 = mr.scaleGrad(g, 0);
    testCase.verifyEqual(g0.amplitude, 0, 'AbsTol', 1e-10);
    testCase.verifyEqual(g0.area, 0, 'AbsTol', 1e-10);
end

%% Test scale arbitrary gradient
function test_scale_arbitrary_grad(testCase)
    sys = mr.opts();
    waveform = [0 10000 20000 10000 0]';
    g = mr.makeArbitraryGrad('x', waveform, sys, 'first', 0, 'last', 0);
    g2 = mr.scaleGrad(g, 0.5);
    testCase.verifyEqual(g2.waveform, 0.5*g.waveform, 'AbsTol', 1e-6);
    testCase.verifyEqual(g2.first, 0.5*g.first, 'AbsTol', 1e-6);
    testCase.verifyEqual(g2.last, 0.5*g.last, 'AbsTol', 1e-6);
end

%% Test id field removed after scaling
function test_id_field_removed(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g.id = 42;  % simulate an assigned library ID
    g2 = mr.scaleGrad(g, 1.5);
    testCase.verifyFalse(isfield(g2, 'id'), 'id field should be removed after scaling');
end

%% Test maxGrad violation with system
function test_maxGrad_violation(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    sys = mr.opts();
    % Scale to exceed maxGrad
    scale_factor = 2 * sys.maxGrad / abs(g.amplitude);
    verifyErrorThrown(testCase, @() mr.scaleGrad(g, scale_factor, sys));
end

%% Test maxSlew violation with system
function test_maxSlew_violation(testCase)
    % Create a gradient with moderate amplitude and short rise time
    g = mr.makeTrapezoid('x', 'amplitude', 1e6, 'riseTime', 1e-4, 'flatTime', 1e-4, 'fallTime', 1e-4);
    sys = mr.opts();
    % Scale to exceed slew rate
    big_scale = 100;
    verifyErrorThrown(testCase, @() mr.scaleGrad(g, big_scale, sys));
end

%% Helper function
function verifyErrorThrown(testCase, funcHandle)
    didError = false;
    try
        funcHandle();
    catch
        didError = true;
    end
    testCase.verifyTrue(didError, 'Expected an error, but none was thrown.');
end
