function tests = testSplitGradientAt
    tests = functiontests(localfunctions);
end

%% Test split trapezoid at flat-top midpoint
function test_split_at_flat_top(testCase)
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3, 'system', sys);
    mid = mr.calcDuration(g) / 2;
    [g1, g2] = mr.splitGradientAt(g, mid, sys);
    testCase.verifyTrue(strcmp(g1.type, 'grad'));
    testCase.verifyTrue(strcmp(g2.type, 'grad'));
    testCase.verifyEqual(g1.channel, 'x');
    testCase.verifyEqual(g2.channel, 'x');
end

%% Test split at rise time
function test_split_at_rise(testCase)
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3, 'system', sys);
    cut_t = g.riseTime;  % at end of ramp-up
    [g1, g2] = mr.splitGradientAt(g, cut_t, sys);
    testCase.verifyTrue(strcmp(g1.type, 'grad'));
    testCase.verifyTrue(strcmp(g2.type, 'grad'));
end

%% Test cut point after gradient throws error
function test_cut_after_end_error(testCase)
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3, 'system', sys);
    total = mr.calcDuration(g);
    verifyErrorThrown(testCase, @() mr.splitGradientAt(g, total + 1e-3, sys));
end

%% Test nargout=1 returns array
function test_single_output(testCase)
    sys = mr.opts();
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3, 'system', sys);
    mid = mr.calcDuration(g) / 2;
    grads = mr.splitGradientAt(g, mid, sys);
    testCase.verifyEqual(length(grads), 2);
    testCase.verifyEqual(mr.calcDuration(grads(1)), grads(2).delay, 'AbsTol', 1e-9); 
    testCase.verifyEqual(mr.calcDuration(grads(1)), mr.calcDuration(grads(2))-grads(2).delay, 'AbsTol', 1e-9); % FIXME: this test will fail if mid is not in grad raster, but 'Duration' above is well divisible by 2
end

%% Test split arbitrary gradient
function test_split_arbitrary(testCase)
    sys = mr.opts();
    w = linspace(0, 20000, 20)';
    w = [w; w(end:-1:1)]; 
    g = mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0);
    mid = mr.calcDuration(g) / 2;
    [g1, g2] = mr.splitGradientAt(g, mid, sys);
    testCase.verifyEqual(g1.channel, 'x');
    testCase.verifyEqual(g2.channel, 'x');
    testCase.verifyEqual([g1.waveform; g2.waveform], w);
    testCase.verifyEqual(g1.last, g2.first, 'AbsTol', 1e-6);
    testCase.verifyEqual(g1.shape_dur+g2.shape_dur, length(w)*sys.gradRasterTime, 'AbsTol', 1e-9);
end

%% Test split arbitrary gradient
function test_split_arbitrary_with_oversampling(testCase)
    sys = mr.opts();
    w = linspace(0, 20000, 21)';
    %w = [w(2:end); w(end); w(end:-1:2)];  % odd count for oversampling, removing the end-zeros because they will be in first/last
    w1 = [w(2:end); w(end-1:-1:2)];  % odd count for oversampling, removing the end-zeros because they will be in first/last
    g = mr.makeArbitraryGrad('x', w1, sys, 'first', 0, 'last', 0, 'oversampling', true);
    mid = mr.calcDuration(g) / 2;
    [g1, g2] = mr.splitGradientAt(g, mid, sys);
    testCase.verifyEqual(g1.channel, 'x');
    testCase.verifyEqual(g1.waveform, w(2:end-1));
    testCase.verifyEqual(g2.waveform, w(end-1:-1:2));
    testCase.verifyEqual(g1.last, g2.first,'AbsTol', 1e-6);
    testCase.verifyEqual(g1.last, w(end),'AbsTol', 1e-6);
    testCase.verifyEqual(g1.shape_dur+g2.shape_dur, (length(w1)+1)*sys.gradRasterTime*0.5, 'AbsTol', 1e-9);
end

%% Test split extended trapezoid
function test_split_extended_trap(testCase)
    sys = mr.opts();
    g = mr.makeExtendedTrapezoid('x', ...
        'Times', [0 2e-4 8e-4 1e-3], ...
        'Amplitudes', [0 50000 50000 0]);
    [g1, g2] = mr.splitGradientAt(g, 5e-4, sys);
    testCase.verifyTrue(strcmp(g1.type, 'grad'));
    testCase.verifyTrue(strcmp(g2.type, 'grad'));
    testCase.verifyEqual(g1.shape_dur+g2.shape_dur, g.shape_dur, 'AbsTol', 1e-9);
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
