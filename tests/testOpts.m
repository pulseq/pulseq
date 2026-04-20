function tests = testOpts
    tests = functiontests(localfunctions);
end

%% Test default opts returns required fields
function test_default_fields(testCase)
    sys = mr.opts();
    required_fields = {'maxGrad', 'maxSlew', 'maxB1', 'riseTime', ...
        'rfDeadTime', 'rfRingdownTime', 'adcDeadTime', ...
        'adcRasterTime', 'rfRasterTime', 'gradRasterTime', ...
        'blockDurationRaster', 'gamma', 'B0'};
    for i = 1:length(required_fields)
        testCase.verifyTrue(isfield(sys, required_fields{i}), ...
            sprintf('Default opts should have field %s', required_fields{i}));
    end
end

%% Test default values
function test_default_values(testCase)
    sys = mr.opts();
    % Default maxGrad = 40 mT/m converted to Hz/m
    expected_maxGrad = mr.convert(40, 'mT/m', 'Hz/m');
    testCase.verifyEqual(sys.maxGrad, expected_maxGrad, 'RelTol', 0.01);
    % Default gradRasterTime = 10 us
    testCase.verifyEqual(sys.gradRasterTime, 10e-6, 'AbsTol', 1e-10);
    % Default rfRasterTime = 1 us
    testCase.verifyEqual(sys.rfRasterTime, 1e-6, 'AbsTol', 1e-10);
end

%% Test custom maxGrad in mT/m
function test_custom_maxGrad(testCase)
    sys = mr.opts('maxGrad', 30, 'gradUnit', 'mT/m');
    expected = mr.convert(30, 'mT/m', 'Hz/m');
    testCase.verifyEqual(sys.maxGrad, expected, 'RelTol', 0.01);
end

%% Test custom maxSlew
function test_custom_maxSlew(testCase)
    sys = mr.opts('maxSlew', 100, 'slewUnit', 'T/m/s');
    expected = mr.convert(100, 'T/m/s', 'Hz/m/s');
    testCase.verifyEqual(sys.maxSlew, expected, 'RelTol', 0.01);
end

%% Test riseTime override
function test_riseTime(testCase)
    sys = mr.opts('riseTime', 250e-6, 'maxGrad', 40, 'gradUnit', 'mT/m');
    testCase.verifyEqual(sys.riseTime, 250e-6);
    % maxSlew should be maxGrad/riseTime
    testCase.verifyEqual(sys.maxSlew, sys.maxGrad / 250e-6, 'RelTol', 0.01);
end

%% Test gamma and B0
function test_gamma_B0(testCase)
    sys = mr.opts();
    testCase.verifyEqual(sys.gamma, 42576000, 'AbsTol', 1);
    testCase.verifyEqual(sys.B0, 1.5, 'AbsTol', 0.01);
end

%% Test custom rfDeadTime
function test_custom_rfDeadTime(testCase)
    sys = mr.opts('rfDeadTime', 100e-6);
    testCase.verifyEqual(sys.rfDeadTime, 100e-6);
end

%% Test repeated calls return same defaults
function test_consistency(testCase)
    sys1 = mr.opts();
    sys2 = mr.opts();
    testCase.verifyEqual(sys1.maxGrad, sys2.maxGrad);
    testCase.verifyEqual(sys1.maxSlew, sys2.maxSlew);
    testCase.verifyEqual(sys1.gradRasterTime, sys2.gradRasterTime);
end
