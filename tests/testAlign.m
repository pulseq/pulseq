function tests = testAlign
    tests = functiontests(localfunctions);
end

%% Test left alignment sets delay to 0
function test_left_align(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    [rf_out, gx_out] = mr.align('left', rf, gx);
    testCase.verifyEqual(rf_out.delay, 0, 'AbsTol', 1e-10);
    testCase.verifyEqual(gx_out.delay, 0, 'AbsTol', 1e-10);
end

%% Test right alignment
function test_right_align(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    [rf_out, gx_out] = mr.align('right', rf, gx);
    % The longer event (gx) should have delay near 0
    % The shorter event (rf) should be right-aligned (delay > 0)
    total_dur_in = mr.calcDuration(rf, gx);
    total_dur_out = mr.calcDuration(rf_out, gx_out);
    testCase.verifyEqual(total_dur_in, total_dur_out, 'AbsTol', 1e-9);
    testCase.verifyEqual(mr.calcDuration(rf_out), total_dur_out, 'AbsTol', 1e-9);
    testCase.verifyEqual(mr.calcDuration(gx_out), total_dur_out, 'AbsTol', 1e-9);
end

%% Test center alignment
function test_center_align(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    rf = mr.makeBlockPulse(pi/2, 'Duration', gx.flatTime);
    [rf_out, gx_out] = mr.align('center', rf, gx);
    % RF should be centered within the block duration
    testCase.verifyEqual(rf_out.delay, gx_out.riseTime, 'AbsTol', 1e-6, 'RF delay should be equal to the ramp-up time of the gradient');
    testCase.verifyEqual(gx_out.delay, 0, 'AbsTol', 1e-6, 'Gradient delay should be zero when RF is centered on it');
end


%% Test with required duration
function test_required_duration(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    required = 3e-3;
    [rf_out, gx_out] = mr.align('right', rf, gx, required);
    testCase.verifyEqual(mr.calcDuration(rf_out), required, 'AbsTol', 1e-6);
    testCase.verifyEqual(mr.calcDuration(gx_out), required, 'AbsTol', 1e-6);
end

%% Test required duration too short throws error
function test_required_duration_too_short(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    verifyErrorThrown(testCase, @() mr.align('left', gx, 0.1e-3));
end

%% Test first parameter must be string
function test_first_param_must_be_string(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    verifyErrorThrown(testCase, @() mr.align(gx, gx));
end

%% Test mixed alignment
function test_mixed_alignment(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 2e-3);
    [rf_out, gx_out] = mr.align('left', rf, 'right', gx);
    testCase.verifyEqual(rf_out.delay, 0, 'AbsTol', 1e-10);
    % gx should be right-aligned
    total_dur = max(mr.calcDuration(rf_out), mr.calcDuration(gx_out));
    testCase.verifyEqual(mr.calcDuration(gx_out), total_dur, 'AbsTol', 1e-6);
    % RF should be left-aligned (delay=0)
    testCase.verifyEqual(rf_out.delay, 0, 'AbsTol', 1e-9);
end

%% Test single output (cell)
function test_single_output(testCase)
    rf = mr.makeBlockPulse(pi/2, 'Duration', 0.5e-3);
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    out = mr.align('left', rf, gx);
    testCase.verifyTrue(iscell(out));
end

%% Helper
function verifyErrorThrown(testCase, funcHandle)
    didError = false;
    try
        funcHandle();
    catch
        didError = true;
    end
    testCase.verifyTrue(didError, 'Expected an error, but none was thrown.');
end
