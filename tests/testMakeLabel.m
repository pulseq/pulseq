function tests = testMakeLabel
    tests = functiontests(localfunctions);
end

%% Test SET label
function test_set_label(testCase)
    lbl = mr.makeLabel('SET', 'SLC', 5);
    testCase.verifyEqual(lbl.type, 'labelset');
    testCase.verifyEqual(lbl.label, 'SLC');
    testCase.verifyEqual(lbl.value, 5);
end

%% Test INC label
function test_inc_label(testCase)
    lbl = mr.makeLabel('INC', 'LIN', 1);
    testCase.verifyEqual(lbl.type, 'labelinc');
    testCase.verifyEqual(lbl.label, 'LIN');
    testCase.verifyEqual(lbl.value, 1);
end

%% Test negative increment
function test_negative_increment(testCase)
    lbl = mr.makeLabel('INC', 'PAR', -1);
    testCase.verifyEqual(lbl.value, -1);
end

%% Test invalid label name throws error
function test_invalid_label_error(testCase)
    testCase.verifyError(@() mr.makeLabel('SET', 'INVALID_LABEL', 0), ...
        'makeLabel:invalidArguments');
end

%% Test invalid type throws error
function test_invalid_type_error(testCase)
    testCase.verifyError(@() mr.makeLabel('RESET', 'SLC', 0), ...
        'makeLabel:invalidArguments');
end

%% Test logical value for flag
function test_logical_flag(testCase)
    lbl = mr.makeLabel('SET', 'REV', true);
    testCase.verifyEqual(lbl.value, true);
end

%% Test missing arguments throws error
function test_missing_args_error(testCase)
    testCase.verifyError(@() mr.makeLabel('SET', 'SLC'), ...
        'makeLabel:invalidArguments');
end

%% Test all standard counters
function test_standard_counters(testCase)
    counters = {'SLC', 'SEG', 'REP', 'AVG', 'SET', 'ECO', 'PHS', 'LIN', 'PAR', 'ACQ'};
    for i = 1:length(counters)
        lbl = mr.makeLabel('SET', counters{i}, 0);
        testCase.verifyEqual(lbl.label, counters{i});
    end
end

%% Test all standard flags
function test_standard_flags(testCase)
    flags = {'NAV', 'REV', 'SMS', 'REF', 'IMA', 'NOISE'};
    for i = 1:length(flags)
        lbl = mr.makeLabel('SET', flags{i}, true);
        testCase.verifyEqual(lbl.label, flags{i});
    end
end
