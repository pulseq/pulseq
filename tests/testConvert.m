function tests = testConvert
    tests = functiontests(localfunctions);
end

%% Test gradient mT/m to Hz/m
function test_gradient_mTm_to_Hzm(testCase)
    gamma = 42.576e6; % Hz/T
    out = mr.convert(1, 'mT/m', 'Hz/m');
    testCase.verifyEqual(out, 1e-3*gamma, 'AbsTol', 1);
end

%% Test slew T/m/s to Hz/m/s
function test_slew_Tms_to_Hzms(testCase)
    gamma = 42.576e6;
    out = mr.convert(1, 'T/m/s', 'Hz/m/s');
    testCase.verifyEqual(out, gamma, 'AbsTol', 1);
end

%% Test B1 T to Hz
function test_b1_T_to_Hz(testCase)
    gamma = 42.576e6;
    out = mr.convert(1, 'T', 'Hz');
    testCase.verifyEqual(out, gamma, 'AbsTol', 1);
end

%% Test B1 uT to Hz
function test_b1_uT_to_Hz(testCase)
    gamma = 42.576e6;
    out = mr.convert(1, 'uT', 'Hz');
    testCase.verifyEqual(out, 1e-6*gamma, 'AbsTol', 0.1);
end

%% Test round-trip gradient
function test_round_trip_gradient(testCase)
    x = 12345;
    out = mr.convert(mr.convert(x, 'Hz/m', 'mT/m'), 'mT/m', 'Hz/m');
    testCase.verifyEqual(out, x, 'AbsTol', 1e-6);
end

%% Test round-trip slew
function test_round_trip_slew(testCase)
    x = 170e6;
    out = mr.convert(mr.convert(x, 'Hz/m/s', 'T/m/s'), 'T/m/s', 'Hz/m/s');
    testCase.verifyEqual(out, x, 'AbsTol', 1e-3);
end

%% Test mT/m/ms to Hz/m/s
function test_slew_mTmms_to_Hzms(testCase)
    gamma = 42.576e6;
    out = mr.convert(1, 'mT/m/ms', 'Hz/m/s');
    testCase.verifyEqual(out, gamma, 'AbsTol', 1);
end

%% Test rad/ms/mm to Hz/m
function test_grad_radmsmm_to_Hzm(testCase)
    out = mr.convert(2*pi*1e-6, 'rad/ms/mm', 'Hz/m');
    testCase.verifyEqual(out, 1, 'AbsTol', 1e-6);
end

%% Test rad/ms/mm/ms to Hz/m/s
function test_slew_radmsmmms_to_Hzms(testCase)
    out = mr.convert(2*pi*1e-9, 'rad/ms/mm/ms', 'Hz/m/s');
    testCase.verifyEqual(out, 1, 'AbsTol', 1e-6);
end

%% Test default toUnit for gradient
function test_default_toUnit_gradient(testCase)
    % Default toUnit for gradients should be 'Hz/m'
    out = mr.convert(1, 'mT/m');
    gamma = 42.576e6;
    testCase.verifyEqual(out, 1e-3*gamma, 'AbsTol', 1);
end

%% Test default toUnit for slew
function test_default_toUnit_slew(testCase)
    out = mr.convert(1, 'T/m/s');
    gamma = 42.576e6;
    testCase.verifyEqual(out, gamma, 'AbsTol', 1);
end

%% Test custom gamma
function test_custom_gamma(testCase)
    custom_gamma = 10000;
    out = mr.convert(1, 'T', 'Hz', 'gamma', custom_gamma);
    testCase.verifyEqual(out, custom_gamma, 'AbsTol', 1e-6);
end

%% Test identity conversion Hz/m to Hz/m
function test_identity_conversion(testCase)
    out = mr.convert(42576, 'Hz/m', 'Hz/m');
    testCase.verifyEqual(out, 42576, 'AbsTol', 1e-10);
end

%% Test mT to Hz
function test_b1_mT_to_Hz(testCase)
    gamma = 42.576e6;
    out = mr.convert(1, 'mT', 'Hz');
    testCase.verifyEqual(out, 1e-3*gamma, 'AbsTol', 1);
end
