function tests = testRestoreAdditionalShapeSamples
    tests = functiontests(localfunctions);
end

%% Test trapezoid-derived shape restores corners
function test_trapezoid_corners(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    % Create a trapezoid and convert to extended-trap-like timing
    g = mr.makeTrapezoid('x', 'Area', 5000, 'Duration', 4e-3);
    % Create a simple shape from the trap (flat portion)
    n = round(g.flatTime / grt);
    waveform = g.amplitude * ones(n, 1);
    tt = ((1:n)' - 0.5) * grt;
    first = 0;
    last = 0;
    [tt_chg, wf_chg] = mr.restoreAdditionalShapeSamples(tt, waveform, first, last, grt, 1);
    testCase.verifyTrue(~isempty(tt_chg));
    testCase.verifyTrue(~isempty(wf_chg));
    testCase.verifyEqual(length(tt_chg), length(wf_chg));
end

%% Test output is valid time-waveform pair
function test_output_valid(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    n = 10;
    waveform = linspace(1000, 10000, n)';
    tt = ((1:n)' - 0.5) * grt;
    first = 500;
    last = 11000;
    warning('OFF', 'mr:restoreShape');
    [tt_chg, wf_chg] = mr.restoreAdditionalShapeSamples(tt, waveform, first, last, grt, 1);
    testCase.verifyTrue(all(isfinite(wf_chg)));
    testCase.verifyTrue(all(isfinite(tt_chg)));
    testCase.verifyTrue(issorted(tt_chg));
end

%% Test constant waveform
function test_constant_waveform(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    n = 20;
    waveform = 5000 * ones(n, 1);
    tt = ((1:n)' - 0.5) * grt;
    first = 5000;
    last = 5000;
    [tt_chg, wf_chg] = mr.restoreAdditionalShapeSamples(tt, waveform, first, last, grt, 1);
    % For a constant waveform, the restored shape should still be constant
    testCase.verifyEqual(length(tt_chg), length(wf_chg));
end
