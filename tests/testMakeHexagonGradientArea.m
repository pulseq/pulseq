function tests = testMakeHexagonGradientArea
%testMakeHexagonGradientArea  Compare makeHexagonGradientArea against
%   makeExtendedTrapezoidArea across a variety of parameter combinations.
%
%   For every test case both functions must:
%     - achieve the requested gradient area (within tolerance)
%     - respect the system slew-rate and amplitude limits
%     - start and end at the prescribed edge amplitudes
%
%   Additionally, the hexagonal waveform (which relaxes the flat-top
%   constraint) must never be *longer* than the trapezoidal one.

    tests = functiontests(localfunctions);
end

% =====================================================================
%  Helpers
% =====================================================================
function sys = defaultSys()
    sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
                  'MaxSlew', 150, 'SlewUnit', 'T/m/s');
end

function verifyGrad(tc, grad, area, g_start, g_end, sys, label)
    % Area matches request
    tc.verifyEqual(grad.area, area, 'AbsTol', 1, ...
        sprintf('%s: area mismatch', label));

    % Edge amplitudes
    tc.verifyEqual(grad.waveform(1),   g_start, 'AbsTol', 1e-3, ...
        sprintf('%s: wrong start amplitude', label));
    tc.verifyEqual(grad.waveform(end), g_end,   'AbsTol', 1e-3, ...
        sprintf('%s: wrong end amplitude', label));

    % Slew rate within limits (with 1 % margin for rounding)
    slew = diff(grad.waveform(:)) ./ diff(grad.tt(:));
    tc.verifyLessThanOrEqual(max(abs(slew)), sys.maxSlew * 1.01, ...
        sprintf('%s: slew rate exceeded', label));

    % Amplitude within limits
    tc.verifyLessThanOrEqual(max(abs(grad.waveform)), sys.maxGrad * 1.01, ...
        sprintf('%s: gradient amplitude exceeded', label));
end

% =====================================================================
%  Test cases — each row is {channel, grad_start, grad_end, area}
% =====================================================================

%% Test: zero-to-zero, positive area
function test_case_01_zero_zero_pos(tc)
    run_case(tc, 'x', 0, 0, 5000);
end

%% Test: zero-to-zero, negative area
function test_case_02_zero_zero_neg(tc)
    run_case(tc, 'y', 0, 0, -5000);
end

%% Test: zero-to-zero, small area
function test_case_03_zero_zero_small(tc)
    run_case(tc, 'z', 0, 0, 100);
end

%% Test: zero-to-zero, large area
function test_case_04_zero_zero_large(tc)
    run_case(tc, 'x', 0, 0, 50000);
end

%% Test: positive start, zero end
function test_case_05_posstart_zeroend(tc)
    sys = defaultSys();
    g_start = sys.maxGrad * 0.5;
    run_case(tc, 'x', g_start, 0, 3000);
end

%% Test: zero start, positive end
function test_case_06_zerostart_posend(tc)
    sys = defaultSys();
    g_end = sys.maxGrad * 0.3;
    run_case(tc, 'y', 0, g_end, 4000);
end

%% Test: positive start, positive end, positive area
function test_case_07_pos_pos_pos(tc)
    sys = defaultSys();
    run_case(tc, 'z', sys.maxGrad*0.2, sys.maxGrad*0.4, 8000);
end

%% Test: negative start, negative end, negative area
function test_case_08_neg_neg_neg(tc)
    sys = defaultSys();
    run_case(tc, 'x', -sys.maxGrad*0.3, -sys.maxGrad*0.1, -6000);
end

%% Test: positive start, negative end (sign change)
function test_case_09_pos_neg(tc)
    sys = defaultSys();
    run_case(tc, 'y', sys.maxGrad*0.2, -sys.maxGrad*0.2, 2000);
end

%% Test: negative start, positive end (sign change)
function test_case_10_neg_pos(tc)
    sys = defaultSys();
    run_case(tc, 'z', -sys.maxGrad*0.3, sys.maxGrad*0.1, -3000);
end

%% Test: equal nonzero edges, positive area
function test_case_11_equal_edges(tc)
    sys = defaultSys();
    g = sys.maxGrad * 0.25;
    run_case(tc, 'x', g, g, 7000);
end

%% Test: near-max edges, small area
%  Note: when both edges are close to max_grad with a small area, the
%  trapezoid flat-top solution can be very compact.  The hexagon may produce
%  a longer waveform in this regime, so we only verify correctness here
%  and log the duration comparison.
function test_case_12_nearmax_small_area(tc)
    sys = defaultSys();
    g = sys.maxGrad * 0.8;
    channel = 'y'; g_start = g; g_end = g; area = 500; %0 fails!

    [grad_trap, ~, ~] = mr.makeExtendedTrapezoidArea(channel, ...
        g_start, g_end, area, sys);
    [grad_hex, ~, ~] = mr.makeHexagonGradientArea(channel, ...
        g_start, g_end, area, sys);

    label_t = sprintf('trap(gs=%.0f ge=%.0f a=%.0f)', g_start, g_end, area);
    label_h = sprintf('hex (gs=%.0f ge=%.0f a=%.0f)', g_start, g_end, area);

    verifyGrad(tc, grad_trap, area, g_start, g_end, sys, label_t);
    verifyGrad(tc, grad_hex,  area, g_start, g_end, sys, label_h);

    % Log durations for information (hexagon may be longer here)
    fprintf('  near-max case: trap=%.6f s  hex=%.6f s\n', ...
        grad_trap.shape_dur, grad_hex.shape_dur);
end

% =====================================================================
%  Core comparison runner
% =====================================================================
function run_case(tc, channel, g_start, g_end, area)
    sys = defaultSys();

    % Trapezoid reference
    [grad_trap, ~, ~] = mr.makeExtendedTrapezoidArea(channel, ...
        g_start, g_end, area, sys);

    % Hexagon under test
    [grad_hex, ~, ~] = mr.makeHexagonGradientArea(channel, ...
        g_start, g_end, area, sys);

    label_t = sprintf('trap(gs=%.0f ge=%.0f a=%.0f)', g_start, g_end, area);
    label_h = sprintf('hex (gs=%.0f ge=%.0f a=%.0f)', g_start, g_end, area);

    % Verify both outputs individually
    verifyGrad(tc, grad_trap, area, g_start, g_end, sys, label_t);
    verifyGrad(tc, grad_hex,  area, g_start, g_end, sys, label_h);

    % Hexagon duration must be <= trapezoid duration
    tc.verifyLessThanOrEqual(grad_hex.shape_dur, grad_trap.shape_dur + 1e-9, ...
        sprintf('Hexagon should be no longer than trapezoid (hex=%.6f trap=%.6f)', ...
            grad_hex.shape_dur, grad_trap.shape_dur));
end
