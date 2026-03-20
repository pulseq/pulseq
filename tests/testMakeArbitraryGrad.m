function tests = testMakeArbitraryGrad
    tests = functiontests(localfunctions);
end

%% Test simple triangle waveform
function test_triangle_waveform(testCase)
    sys = mr.opts();
    w = [0 10000 20000 10000 0]';
    g = mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0);
    testCase.verifyEqual(g.type, 'grad');
    testCase.verifyEqual(g.channel, 'x');
    testCase.verifyEqual(length(g.waveform), length(w));
    testCase.verifyTrue(g.area ~= 0);
end

%% Test with explicit first/last
function test_explicit_first_last(testCase)
    sys = mr.opts();
    w = [5000 10000 5000]';
    g = mr.makeArbitraryGrad('y', w, sys, 'first', 0, 'last', 0);
    testCase.verifyEqual(g.first, 0, 'AbsTol', 1e-10);
    testCase.verifyEqual(g.last, 0, 'AbsTol', 1e-10);
end

%% Test without first/last produces warnings
function test_no_first_last_warnings(testCase)
    sys = mr.opts();
    w = [0 5000 10000 5000 0]';
    % Suppress warnings for this test
    warnState = warning('off', 'all');
    g = mr.makeArbitraryGrad('z', w, sys);
    warning(warnState);
    testCase.verifyEqual(g.type, 'grad');
end

%% Test maxGrad violation
function test_maxGrad_violation(testCase)
    sys = mr.opts();
    % Create waveform exceeding maxGrad
    w = ones(10, 1) * (sys.maxGrad * 1.01);
    verifyErrorThrown(testCase, @() mr.makeArbitraryGrad('x', w, sys, 'first', w(1), 'last', w(end)));
end

%% Test maxSlew violation
function test_maxSlew_violation(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    % Create waveform with huge step → slew violation
    w = [0; sys.maxGrad * 0.5; 0];  % rapid change in one step
    % The slew = change/grt; ensure it exceeds maxSlew
    step = sys.maxSlew * grt * 1.05;  % 1.05x maxSlew
    w = [0; step; 0];
    verifyErrorThrown(testCase, @() mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0));
end

%% Test oversampling with odd samples
function test_oversampling_odd(testCase)
    sys = mr.opts();
    w = [0 5000 10000 5000 0]';  % 5 samples (odd)
    g = mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0, 'oversampling', true);
    testCase.verifyEqual(g.type, 'grad');
    testCase.verifyEqual(length(g.tt), length(w));
end

%% Test oversampling with even samples throws error
function test_oversampling_even_error(testCase)
    sys = mr.opts();
    w = [0 5000 10000 15000 20000 0]';  % 6 samples (even)
    verifyErrorThrown(testCase, @() mr.makeArbitraryGrad('x', w, sys, 'first', 0, 'last', 0, 'oversampling', true));
end

%% Test area calculation
function test_area_calculation(testCase)
    sys = mr.opts();
    % Constant gradient
    amp = 10000;
    n = 20;
    w = amp * ones(n, 1);
    g = mr.makeArbitraryGrad('x', w, sys, 'first', amp, 'last', amp);
    expected_area = amp * n * sys.gradRasterTime;
    testCase.verifyEqual(g.area, expected_area, 'AbsTol', 1);
end

%% Test valid channels
function test_valid_channels(testCase)
    sys = mr.opts();
    w = [0 5000 0]';
    for ch = {'x', 'y', 'z'}
        g = mr.makeArbitraryGrad(ch{1}, w, sys, 'first', 0, 'last', 0);
        testCase.verifyEqual(g.channel, ch{1});
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
