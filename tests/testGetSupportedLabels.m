function tests = testGetSupportedLabels
    tests = functiontests(localfunctions);
end

%% Setup: reset global state
function setup(~)
    mr.aux.globalVars('set', 'SupportedLabels', []);
end

%% Teardown: reset global state
function teardown(~)
    mr.aux.globalVars('set', 'SupportedLabels', []);
end

%% Test default labels contain expected entries
function test_default_labels(testCase)
    labels = mr.getSupportedLabels();
    testCase.verifyTrue(iscell(labels));
    expected = {'SLC', 'SEG', 'REP', 'AVG', 'SET', 'ECO', 'PHS', 'LIN', 'PAR', 'ACQ', 'TRID'};
    for i = 1:length(expected)
        testCase.verifyTrue(any(strcmp(labels, expected{i})), ...
            sprintf('Label %s should be in default set', expected{i}));
    end
end

%% Test default labels have at least 23 entries
function test_default_labels_count(testCase)
    labels = mr.getSupportedLabels();
    testCase.verifyGreaterThanOrEqual(length(labels), 23);
end

%% Test multiple calls return same result
function test_consistency(testCase)
    labels1 = mr.getSupportedLabels();
    labels2 = mr.getSupportedLabels();
    testCase.verifyEqual(labels1, labels2);
end

%% Test flags are present
function test_flags_present(testCase)
    labels = mr.getSupportedLabels();
    flags = {'NAV', 'REV', 'SMS', 'REF', 'IMA', 'NOISE'};
    for i = 1:length(flags)
        testCase.verifyTrue(any(strcmp(labels, flags{i})), ...
            sprintf('Flag %s should be in default set', flags{i}));
    end
end

%% Test control labels present
function test_control_labels_present(testCase)
    labels = mr.getSupportedLabels();
    ctrl = {'PMC', 'NOROT', 'NOPOS', 'NOSCL', 'ONCE'};
    for i = 1:length(ctrl)
        testCase.verifyTrue(any(strcmp(labels, ctrl{i})), ...
            sprintf('Control label %s should be in default set', ctrl{i}));
    end
end
