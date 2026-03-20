function tests = testCalcRfCenter
    tests = functiontests(localfunctions);
end

%% Test block pulse center at midpoint
function test_block_pulse_center(testCase)
    dur = 1e-3;
    rf = mr.makeBlockPulse(pi/2, 'Duration', dur);
    [tc, ic, fi] = mr.calcRfCenter(rf);
    testCase.verifyEqual(tc, dur/2, 'AbsTol', 1e-6, ...
        'Block pulse center should be at midpoint');
    testCase.verifyTrue(ic >= 1, 'Integer index should be positive');
    testCase.verifyTrue(abs(fi) <= 0.5, 'Fractional index should be in [-0.5, 0.5]');
end

%% Test sinc pulse center near peak
function test_sinc_pulse_center(testCase)
    dur = 4e-3;
    rf = mr.makeSincPulse(pi/2, 'Duration', dur, 'timeBwProduct', 4);
    [tc, ~, ~] = mr.calcRfCenter(rf);
    % Sinc peak is typically at center of duration
    testCase.verifyEqual(tc, dur/2, 'AbsTol', 1e-6, ...
        'Sinc pulse center should be near duration/2');
end

%% Test pulse with explicit center field
function test_explicit_center(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 1e-3);
    rf.center = 0.3e-3;
    [tc, ~, ~] = mr.calcRfCenter(rf);
    testCase.verifyEqual(tc, 0.3e-3, 'AbsTol', 1e-9, ...
        'Should return the explicit center value');
end

%% Test three outputs
function test_three_outputs(testCase)
    rf = mr.makeBlockPulse(pi/4, 'Duration', 0.5e-3);
    [tc, ic, fi] = mr.calcRfCenter(rf);
    testCase.verifyTrue(isscalar(tc));
    testCase.verifyTrue(isscalar(ic));
    testCase.verifyTrue(isscalar(fi));
end

%% Test gauss pulse center
function test_gauss_pulse_center(testCase)
    dur = 3e-3;
    rf = mr.makeGaussPulse(pi/2, 'Duration', dur);
    [tc, ~, ~] = mr.calcRfCenter(rf);
    % Gaussian pulse should peak in the middle
    testCase.verifyEqual(tc, dur/2, 'AbsTol', 1e-6, ...
        'Gaussian pulse center should be at duration/2');
end
