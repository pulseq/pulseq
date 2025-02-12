function tests = testMakeTrapezoid
    tests = functiontests(localfunctions);
end

function test_make_trapezoid(testCase)
    % Test for missing required parameters
    try
        mr.makeTrapezoid('x');
        error('Expected error for missing parameters not thrown');
    catch e
        expectedErrMsg = "Must supply either 'area', 'flatArea' or 'amplitude', and only one of the three may be specified";
        assert(strcmp(e.message, expectedErrMsg), 'Error message for missing parameters not as expected');
    end

    % Test for flatTime without corresponding parameters
    try
        mr.makeTrapezoid('x', 'flatTime', 10, 'area', 10);
        error('Expected error for flatTime without required parameters not thrown');
    catch e
        expectedErrMsg = "When 'flatTime' is provided either 'flatArea' or 'amplitude' must be provided as well; you may consider providing 'duration', 'area' and optionally ramp times instead.";
        assert(strcmp(e.message, expectedErrMsg), 'Error message for flatTime mismatch not as expected');
    end

    % Test for area too large
    try
        mr.makeTrapezoid('x', 'area', 1e6, 'duration', 1e-6);
        error('Expected area too large error not thrown');
    catch e
        expectedErrMsg = "Requested area is too large for this gradient. Minimum required duration for this area (accounting for the gradient raster time) is ";
        assert(contains(e.message, expectedErrMsg), 'Error message for area too large not as expected');
    end

    % Test for no area and no duration
    try
        mr.makeTrapezoid('x', 'amplitude', 1);
        error('Expected error for missing area or duration not thrown');
    catch e
        assert(strcmp(e.message, 'Must supply area or duration'), 'Error message for missing area/duration not as expected');
    end

    % Test for amplitude too large
    try
        mr.makeTrapezoid('x', 'amplitude', 1e10, 'duration', 1);
        error('Expected amplitude too large error not thrown');
    catch e
        assert(contains(e.message, 'Amplitude violation'), 'Error message for amplitude too large not as expected');
    end

    % Test for duration too short
    try
        mr.makeTrapezoid('x', 'area', 1, 'duration', 0.1, 'riseTime', 0.1);
        error('Expected duration too short error not thrown');
    catch e
        assert(contains(e.message, 'Requested area is too large for this gradient duration. Probably amplitude is violated'), 'Error message for short duration not as expected');
    end

    % Test for minimum input cases
    opts = mr.opts;
    trap = mr.makeTrapezoid('x', 'amplitude', 1, 'duration', 1);
    compare_trap_out(trap, 1, 1e-5, 1 - 2e-5, 1e-5);

    % flatTime + amplitude
    trap = mr.makeTrapezoid('x', 'flatTime', 1, 'amplitude', 1);
    compare_trap_out(trap, 1, 1e-5, 1, 1e-5);

    % flatArea + flatTime
    trap = mr.makeTrapezoid('x', 'flatTime', 1, 'flatArea', 1);
    compare_trap_out(trap, 1, 1e-5, 1, 1e-5);

    % area specified
    trap = mr.makeTrapezoid('x', 'area', 1);
    compare_trap_out(trap, 5e4, 2e-5, 0, 2e-5);
    
    % area + duration + riseTime
    trap = mr.makeTrapezoid('x', 'area', 1, 'duration', 1, 'riseTime', 0.01);
    compare_trap_out(trap, 1 / 0.99, 0.01, 0.98, 0.01);

end

% Helper function for comparing trapezoid output
function compare_trap_out(trap, amplitude, riseTime, flatTime, fallTime)
    assert(abs(trap.amplitude - amplitude) < 1e-10, 'Amplitude mismatch');
    assert(abs(trap.riseTime - riseTime) < 1e-10, 'Rise time mismatch');
    assert(abs(trap.flatTime - flatTime) < 1e-10, 'Flat time mismatch');
    assert(abs(trap.fallTime - fallTime) < 1e-10, 'Fall time mismatch');
end

% Helper function for rounding to raster time
function y = round2raster(x, raster)
    if nargin < 2
        raster = 1e-5;
    end
    y = round(x / raster) * raster;
end
