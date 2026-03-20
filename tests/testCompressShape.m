function tests = testCompressShape
    tests = functiontests(localfunctions);
end

%% Test constant waveform compresses well
function test_constant_waveform(testCase)
    w = ones(1, 100);
    s = mr.compressShape(w);
    testCase.verifyEqual(s.num_samples, 100);
    testCase.verifyTrue(length(s.data) < length(w), ...
        'Constant waveform should compress significantly');
end

%% Test linear ramp compresses well
function test_linear_ramp(testCase)
    w = linspace(0, 1, 100);
    s = mr.compressShape(w);
    testCase.verifyEqual(s.num_samples, 100);
    testCase.verifyTrue(length(s.data) < length(w), ...
        'Linear ramp should compress (constant derivative)');
end

%% Test random waveform may not compress
function test_random_waveform(testCase)
    rng(42);
    w = randn(1, 100);
    s = mr.compressShape(w);
    testCase.verifyEqual(s.num_samples, 100);
    % Random data: compressed form stored only if shorter
    testCase.verifyEqual(length(s.data), s.num_samples, ...
        'Random waveform should not benefit from compression');
end

%% Test short waveform stored as-is by default
function test_short_waveform_no_compression(testCase)
    w = [1 2 3 4];
    s = mr.compressShape(w);
    testCase.verifyEqual(s.num_samples, 4);
    testCase.verifyEqual(s.data, w, 'AbsTol', 1e-10);
end

%% Test short waveform with forceCompression
function test_short_waveform_force_compression(testCase)
    w = [1 1 1 1];
    s = mr.compressShape(w, true);
    testCase.verifyEqual(s.num_samples, 4);
    % forced compression: data should differ from raw
    testCase.verifyTrue(length(s.data) <= length(w));
end

%% Test non-finite sample throws error
function test_inf_sample_error(testCase)
    w = [1 2 Inf 4];
    testCase.verifyError(@() mr.compressShape(w), ?MException);
end

function test_nan_sample_error(testCase)
    w = [1 NaN 3 4 5];
    testCase.verifyError(@() mr.compressShape(w), ?MException);
end

%% Test num_samples field always correct
function test_num_samples_field(testCase)
    lengths = [1, 2, 5, 50, 200];
    for n = lengths
        w = sin(linspace(0, 2*pi, n));
        s = mr.compressShape(w);
        testCase.verifyEqual(s.num_samples, n);
    end
end

%% Test round-trip with decompressShape
function test_round_trip(testCase)
    w = sin(linspace(0, 4*pi, 200));
    s = mr.compressShape(w);
    w2 = mr.decompressShape(s);
    testCase.verifyEqual(w2, w(:), 'AbsTol', 1e-6);
end
