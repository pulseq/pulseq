function tests = testPts2waveform
    tests = functiontests(localfunctions);
end

%% Test two-point linear ramp
function test_two_point_ramp(testCase)
    gradRasterTime = 10e-6;
    times = [0 1e-3];
    amplitudes = [0 1000];
    w = mr.pts2waveform(times, amplitudes, gradRasterTime);
    % Output should have approx 100 samples (1ms / 10us)
    expected_len = round((max(times) - min(times)) / gradRasterTime) - 1;
    testCase.verifyEqual(length(w), expected_len, 'AbsTol', 1);
    % First sample at raster center: 0.5*gradRasterTime → amplitude ≈ 5
    testCase.verifyEqual(w(1), 1000 * 0.5 * gradRasterTime / 1e-3, 'AbsTol', 1e-6);
end

%% Test constant amplitude
function test_constant_amplitude(testCase)
    gradRasterTime = 10e-6;
    times = [0 5e-4];
    amplitudes = [42000 42000];
    w = mr.pts2waveform(times, amplitudes, gradRasterTime);
    testCase.verifyEqual(w, 42000 * ones(size(w)), 'AbsTol', 1e-6);
end

%% Test multi-point piecewise linear
function test_multi_point(testCase)
    gradRasterTime = 10e-6;
    times = [0 1e-4 2e-4 3e-4];
    amplitudes = [0 10000 10000 0];
    w = mr.pts2waveform(times, amplitudes, gradRasterTime);
    % Should be a valid waveform
    testCase.verifyTrue(~isempty(w));
    testCase.verifyTrue(all(isfinite(w)));
end

%% Test output is column or row vector
function test_output_vector(testCase)
    gradRasterTime = 10e-6;
    times = [0 2e-4];
    amplitudes = [0 5000];
    w = mr.pts2waveform(times, amplitudes, gradRasterTime);
    % Should be a numeric vector
    testCase.verifyTrue(isnumeric(w));
    testCase.verifyTrue(isvector(w));
end
