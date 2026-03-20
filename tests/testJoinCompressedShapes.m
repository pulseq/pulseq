function tests = testJoinCompressedShapes
    tests = functiontests(localfunctions);
end

%% Test num_samples is sum of inputs
function test_num_samples_sum(testCase)
    w1 = [0 1 2 1 0]';
    w2 = [0 3 3 0]';
    s1 = mr.compressShape(w1, true);
    s2 = mr.compressShape(w2, true);
    s = mr.joinCompressedShapes(s1, s2);
    testCase.verifyEqual(s.num_samples, s1.num_samples + s2.num_samples);
end

%% Test equivalence to decompress-concat-recompress
function test_join_equivalence(testCase)
    w1 = [0 1 2 1 0]';
    w2 = [0 3 3 0]';
    s1 = mr.compressShape(w1, true);
    s2 = mr.compressShape(w2, true);
    s = mr.joinCompressedShapes(s1, s2);
    w_joined = mr.decompressShape(s);
    w_expected = [w1; w2];
    testCase.verifyEqual(w_joined, w_expected, 'AbsTol', 1e-6);
end

%% Test joining shapes with zero runs at boundary
function test_join_zero_runs(testCase)
    % Create shapes with runs of zeros at the boundary
    w1 = [0 1 0 0 0 0]';  % ends with zero run
    w2 = [0 0 0 0 2 0]';  % starts with zero run
    s1 = mr.compressShape(w1, true);
    s2 = mr.compressShape(w2, true);
    s = mr.joinCompressedShapes(s1, s2);
    w_joined = mr.decompressShape(s);
    testCase.verifyEqual(w_joined, [w1; w2], 'AbsTol', 1e-6);
    testCase.verifyEqual(s.num_samples, 12);
end

%% Test simple concatenation case
function test_simple_concat(testCase)
    w1 = [0 1 2 0]';
    w2 = [0 5 6 0]';
    s1 = mr.compressShape(w1, true);
    s2 = mr.compressShape(w2, true);
    s = mr.joinCompressedShapes(s1, s2);
    w_joined = mr.decompressShape(s);
    testCase.verifyEqual(w_joined, [w1; w2], 'AbsTol', 1e-6);
end

%% Test joining longer shapes
function test_join_long_shapes(testCase)
    w1 = [zeros(1,50) linspace(0,1,50) zeros(1,50)]';
    w2 = [zeros(1,30) ones(1,40) zeros(1,30)]';
    s1 = mr.compressShape(w1);
    s2 = mr.compressShape(w2);
    s = mr.joinCompressedShapes(s1, s2);
    w_joined = mr.decompressShape(s);
    testCase.verifyEqual(w_joined, [w1; w2], 'AbsTol', 1e-6);
end
