function tests = testMakeArbitraryRf
    tests = functiontests(localfunctions);
end

%% Test rectangular signal with pi/2 flip
function test_rectangular_signal(testCase)
    sys = mr.opts();
    n = 100;
    signal = ones(1, n);
    flip = pi/6;
    rf = mr.makeArbitraryRf(signal, flip, sys);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(length(rf.signal), n);
    % Verify flip angle: sum(signal)*dwell * 2*pi should equal flip
    actual_flip = abs(sum(rf.signal .* sys.rfRasterTime)) * 2 * pi;
    testCase.verifyEqual(actual_flip, flip, 'RelTol', 0.01);
end

%% Test complex signal
function test_complex_signal(testCase)
    sys = mr.opts();
    n = 500;
    signal = ones(1, n) .* exp(1i * linspace(0, pi, n));
    rf = mr.makeArbitraryRf(signal, pi/4, sys);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(length(rf.signal), n);
end

%% Test with explicit center
function test_explicit_center(testCase)
    sys = mr.opts();
    signal = ones(1, 100);
    rf = mr.makeArbitraryRf(signal, pi/6, sys, 'center', 30e-6);
    testCase.verifyEqual(rf.center, 30e-6, 'AbsTol', 1e-9);
end

%% Test bandwidth + sliceThickness returns gradient
function test_with_slice_selection(testCase)
    sys = mr.opts();
    signal = ones(1, 100);
    [rf, gz] = mr.makeArbitraryRf(signal, pi/6, sys, ...
        'bandwidth', 1000, 'sliceThickness', 5e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(gz.type, 'trap');
    testCase.verifyEqual(gz.channel, 'z');
end

%% Test both bandwidth and timeBwProduct error
function test_bw_and_tbw_error(testCase)
    sys = mr.opts();
    signal = ones(1, 100);
    verifyErrorThrown(testCase, @() mr.makeArbitraryRf(signal, pi/6, sys, ...
        'bandwidth', 1000, 'timeBwProduct', 4, 'sliceThickness', 5e-3));
end

%% Test shape_dur matches expected
function test_shape_dur(testCase)
    sys = mr.opts();
    n = 200;
    signal = ones(1, n);
    rf = mr.makeArbitraryRf(signal, pi/6, sys);
    expected_dur = n * sys.rfRasterTime;
    testCase.verifyEqual(rf.shape_dur, expected_dur, 'AbsTol', 1e-10);
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
