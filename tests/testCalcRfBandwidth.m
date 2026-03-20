function tests = testCalcRfBandwidth
    tests = functiontests(localfunctions);
end

%% Test block pulse bandwidth approx 1/duration
function test_block_pulse_bw(testCase)
    dur = 1e-3;
    rf = mr.makeBlockPulse(pi/2, 'Duration', dur);
    [bw, fc] = mr.calcRfBandwidth(rf);
    % Block pulse has BW ≈ 1/duration at 0.5 cutoff (sinc-like spectrum)
    expected_bw = 1/dur;
    testCase.verifyEqual(bw, expected_bw, 'RelTol', 0.3, ...
        'Block pulse BW should be approximately 1/duration');
    % Center frequency should be near zero for on-resonance
    testCase.verifyEqual(fc, 0, 'AbsTol', 50);
end

%% Test sinc pulse bandwidth approx TBW/duration
function test_sinc_pulse_bw(testCase)
    dur = 4e-3;
    tbw = 4;
    rf = mr.makeSincPulse(pi/2, 'Duration', dur, 'timeBwProduct', tbw);
    [bw, ~] = mr.calcRfBandwidth(rf);
    expected_bw = tbw / dur;
    testCase.verifyEqual(bw, expected_bw, 'RelTol', 0.15, ...
        'Sinc pulse BW should be approximately TBW/duration');
end

%% Test custom cutoff
function test_custom_cutoff(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 1e-3);
    [bw_half, ~] = mr.calcRfBandwidth(rf, 0.5);
    [bw_quarter, ~] = mr.calcRfBandwidth(rf, 0.25);
    testCase.verifyTrue(bw_quarter > bw_half, ...
        'Lower cutoff should yield wider bandwidth');
end

%% Test all 6 outputs returned
function test_six_outputs(testCase)
    rf = mr.makeBlockPulse(pi/6, 'Duration', 0.5e-3);
    [bw, fc, spectrum, f, rfs, t] = mr.calcRfBandwidth(rf);
    testCase.verifyTrue(isscalar(bw));
    testCase.verifyTrue(isscalar(fc));
    testCase.verifyTrue(isvector(spectrum));
    testCase.verifyTrue(isvector(f));
    testCase.verifyTrue(isvector(rfs));
    testCase.verifyTrue(isvector(t));
    testCase.verifyEqual(length(spectrum), length(f));
    testCase.verifyEqual(length(rfs), length(t));
end
