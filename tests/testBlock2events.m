function tests = testBlock2events
    tests = functiontests(localfunctions);
end

%% Test block struct with multiple events
function test_block_struct(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 1e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    block = struct('rf', rf, 'gx', gx, 'gy', [], 'gz', []);
    c = mr.block2events(block);
    testCase.verifyTrue(iscell(c));
    % empty fields removed → only rf and gx remain
    testCase.verifyEqual(length(c), 2);
end

%% Test cell array passthrough
function test_cell_passthrough(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    adc = mr.makeAdc(128, 'Duration', 1e-3);
    c = mr.block2events({gx, adc});
    testCase.verifyTrue(iscell(c));
    testCase.verifyEqual(length(c), 2);
end

%% Test nested cell unwrapping
function test_nested_cell_unwrap(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    nested = {{gx}};
    c = mr.block2events(nested);
    % Should unwrap the nested 1x1 cell wrappers
    testCase.verifyTrue(isstruct(c) || (iscell(c) && length(c) == 1));
end

%% Test single event in block
function test_single_event_block(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    block = struct('rf', [], 'gx', gx, 'gy', [], 'gz', []);
    c = mr.block2events(block);
    testCase.verifyTrue(iscell(c));
    testCase.verifyEqual(length(c), 1);
end
