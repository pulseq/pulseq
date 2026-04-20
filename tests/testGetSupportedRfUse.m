function tests = testGetSupportedRfUse
    tests = functiontests(localfunctions);
end

%% Test returns cell array
function test_returns_cell(testCase)
    uses = mr.getSupportedRfUse();
    testCase.verifyTrue(iscell(uses));
end

%% Test contains expected entries
function test_expected_entries(testCase)
    uses = mr.getSupportedRfUse();
    expected = {'excitation', 'refocusing', 'inversion', 'saturation', 'preparation'};
    for i = 1:length(expected)
        testCase.verifyTrue(any(strcmp(uses, expected{i})), ...
            sprintf('Should contain %s', expected{i}));
    end
end

%% Test two outputs
function test_two_outputs(testCase)
    [uses, short] = mr.getSupportedRfUse();
    testCase.verifyTrue(iscell(uses));
    testCase.verifyTrue(ischar(short));
    testCase.verifyEqual(length(short), length(uses));
end

%% Test short names are first letters
function test_short_names(testCase)
    [uses, short] = mr.getSupportedRfUse();
    for i = 1:length(uses)
        testCase.verifyEqual(short(i), uses{i}(1));
    end
end

%% Test underfined and other included
function test_undefined_and_other(testCase)
    uses = mr.getSupportedRfUse();
    testCase.verifyTrue(any(strcmp(uses, 'undefined')));
    testCase.verifyTrue(any(strcmp(uses, 'other')));
end
