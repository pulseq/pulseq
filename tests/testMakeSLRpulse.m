function tests = testMakeSLRpulse
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

%% Test basic SLR pulse creation
function test_basic_slr(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping SLR pulse tests');
    rf = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, 'timeBwProduct', 4);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test SLR with slice thickness returns gradient
function test_slr_with_slice(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping SLR pulse tests');
    [rf, gz] = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, ...
        'timeBwProduct', 4, 'sliceThickness', 5e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(gz.type, 'trap');
end

%% Test SLR excitation pulse
function test_slr_excitation(testCase)
    assumeTrue(testCase, testCase.TestData.pythonAvailable, ...
        'Python/sigpy not available, skipping SLR pulse tests');
    rf = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, ...
        'timeBwProduct', 4, 'use', 'excitation');
    testCase.verifyEqual(rf.use, 'excitation');
end
