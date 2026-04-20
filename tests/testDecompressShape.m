function tests = testDecompressShape
    tests = functiontests(localfunctions);
end

%% Test round-trip for constant waveform
function test_round_trip_constant(testCase)
    w = ones(100, 1);
    s = mr.compressShape(w);
    w2 = mr.decompressShape(s);
    testCase.verifyEqual(w2, w, 'AbsTol', 1e-6);
end

%% Test round-trip for linear ramp
function test_round_trip_linear(testCase)
    w = linspace(0, 1, 100)';
    s = mr.compressShape(w);
    w2 = mr.decompressShape(s);
    testCase.verifyEqual(w2, w, 'AbsTol', 1e-6);
end

%% Test round-trip for sinusoid
function test_round_trip_sinusoid(testCase)
    w = sin(linspace(0, 4*pi, 300))';
    s = mr.compressShape(w);
    w2 = mr.decompressShape(s);
    testCase.verifyEqual(w2, w, 'AbsTol', 1e-6);
end

%% Test round-trip for random waveform
function test_round_trip_random(testCase)
    rng(123);
    w = randn(100, 1);
    s = mr.compressShape(w);
    w2 = mr.decompressShape(s);
    testCase.verifyEqual(w2, w, 'AbsTol', 1e-6);
end

%% Test uncompressed passthrough
function test_uncompressed_passthrough(testCase)
    % If num_samples == length(data), it's uncompressed → passthrough
    s.num_samples = 5;
    s.data = [1 2 3 4 5];
    w = mr.decompressShape(s);
    testCase.verifyEqual(w, [1; 2; 3; 4; 5]);
end

%% Test forceDecompression on uncompressed data
function test_force_decompression(testCase)
    w_orig = linspace(0, 1, 50)';
    % Now make an "uncompressed" shape that is actually compressed
    s_forced = mr.compressShape(w_orig, true);
    w2 = mr.decompressShape(s_forced, true);
    testCase.verifyEqual(w2, w_orig, 'AbsTol', 1e-6);
end

%% Test output length matches num_samples
function test_output_length(testCase)
    for n = [10, 50, 200]
        w = sin(linspace(0, 2*pi, n))';
        s = mr.compressShape(w);
        w2 = mr.decompressShape(s);
        testCase.verifyEqual(length(w2), s.num_samples);
    end
end
