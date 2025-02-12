function tests = testMakeBlockPulse
    tests = functiontests(localfunctions);
end

function test_invalid_use_error(testCase)
    try
        mr.makeBlockPulse(pi, 'duration', 1e-3, 'use', 'foo');
    catch ME
        assert(contains(ME.message, "value of 'use' is invalid"));
    end
end

function test_bandwidth_and_duration_error(testCase)
    try
        pulse = mr.makeBlockPulse(pi);
    catch ME
        % Check for the user warning
        assert(strcmp(ME.message, 'Either bandwidth or duration must be defined'));
    end
end

function test_generation_methods(testCase)
    % Test minimum input cases

    % Case 1: With duration
    pulse = mr.makeBlockPulse(pi, 'duration', 1e-3);
    assert(isstruct(pulse));
    assert(pulse.shape_dur == 1e-3);

    % Case 2: With bandwidth
    pulse = mr.makeBlockPulse(pi, 'bandwidth', 1e3);
    assert(isstruct(pulse));
    assert(pulse.shape_dur == 1 / (4 * 1e3));

    % Case 3: With bandwidth and time_bw_product
    pulse = mr.makeBlockPulse(pi, 'bandwidth', 1e3, 'timeBwProduct', 5);
    assert(isstruct(pulse));
    assert(pulse.shape_dur == 5 / 1e3);
end

function test_amp_calculation(testCase)
    % A 1 ms 180-degree pulse requires 500 Hz gamma B1
    pulse = mr.makeBlockPulse(pi , 'duration', 1e-3);
    assert(abs(pulse.signal(end) - 500) < 1e-3);

    % A 1 ms 90-degree pulse requires 250 Hz gamma B1
    pulse = mr.makeBlockPulse(0.5 * pi, 'duration', 1e-3);
    assert(abs(pulse.signal(end) - 250) < 1e-3);

    % A 2 ms 90-degree pulse requires 125 Hz gamma B1
    pulse = mr.makeBlockPulse(0.5 * pi, 'duration', 2e-3);
    assert(abs(pulse.signal(end) - 125) < 1e-3);
end
