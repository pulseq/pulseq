function tests = testMakeAdiabaticPulse
    tests = functiontests(localfunctions);
end

%% Setup: check Python/sigpy availability
function setup(testCase)
    [status, ~] = system('python3 -c "import sigpy" 2>/dev/null');
    if status ~= 0
        [status, ~] = system('python -c "import sigpy" 2>/dev/null');
    end
    testCase.TestData.pythonAvailable = (status == 0);
end

%% Test create hypsec pulse
function test_hypsec(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    rf = mr.makeAdiabaticPulse('hypsec');
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test create wurst pulse
function test_wurst(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    rf = mr.makeAdiabaticPulse('wurst','duration', 40e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test hypsec with custom duration
function test_hypsec_duration(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    dur = 15e-3;
    rf = mr.makeAdiabaticPulse('hypsec', 'Duration', dur);
    testCase.verifyEqual(rf.shape_dur, dur, 'RelTol', 0.01);
end

%% Test with slice thickness returns gradient
function test_with_slice_thickness(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    [rf, gz] = mr.makeAdiabaticPulse('hypsec', 'sliceThickness', 5e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(gz.type, 'trap');
    testCase.verifyEqual(gz.channel, 'z');
end
