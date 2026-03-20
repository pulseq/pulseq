function tests = testSplitGradient
    tests = functiontests(localfunctions);
end

%% Test split standard trapezoid into 3 parts
function test_split_standard_trap(testCase)
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3);
    grads = mr.splitGradient(g);
    testCase.verifyEqual(length(grads), 3, 'Should produce 3 gradient parts');
    % All should be on channel x
    for i = 1:3
        testCase.verifyEqual(grads(i).channel, 'x');
    end
end

%% Test durations of split parts
function test_split_durations(testCase)
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3, 'system', sys);
    grads = mr.splitGradient(g, sys);
    testCase.verifyEqual(length(grads), 3);
    % The sum of the shape durations should equal the original duration
    total_dur = mr.calcDuration(g);
    % When reconstructed with addGradients the total duration should match
    g_recon = mr.addGradients({grads(1), grads(2), grads(3)}, sys);
    testCase.verifyEqual(mr.calcDuration(g_recon), total_dur, 'AbsTol', 1e-6);
    testCase.verifyEqual(g_recon.area, g.area, 'RelTol', 0.001);
end

%% Test triangular gradient (flatTime=0)
function test_split_triangular(testCase)
    % Create a triangular gradient by specifying small area and long rise time
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'amplitude', 100000, ...
        'riseTime', 5e-4, 'flatTime', 0, 'fallTime', 5e-4, 'system', sys);
    grads = mr.splitGradient(g, sys);
    testCase.verifyEqual(length(grads), 2);
end

%% Test arbitrary gradient throws error
function test_arbitrary_grad_error(testCase)
    sys = mr.opts();
    w = [0 10000 20000 10000 0]';
    g = mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0, 'system', sys);
    verifyErrorThrown(testCase, @() mr.splitGradient(g));
end

%% Test split preserves channel
function test_channel_preserved(testCase)
    g = mr.makeTrapezoid('y', 'Area', 3000, 'Duration', 3e-3);
    grads = mr.splitGradient(g);
    for i = 1:3
        testCase.verifyEqual(grads(i).channel, 'y');
    end
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
