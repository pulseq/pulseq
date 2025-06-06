function tests = testMakeGaussPulse
    tests = functiontests(localfunctions);
end

function test_make_gauss_pulse(testCase)
    % Test the 'use' parameter for valid and invalid inputs
    try
        mr.makeGaussPulse(1, 'use', 'invalid');
        error('Expected error for invalid use parameter not thrown');
    catch ME
        assert(contains(ME.message, "value of 'use' is invalid"));
    end

    % Loop through all supported rf uses and check if make_gauss_pulse works
    supported_rf_uses = mr.getSupportedRfUse(); % Assuming this function returns a list of supported uses
    
    for i = 1:numel(supported_rf_uses)
        use = supported_rf_uses{i};
        pulse = mr.makeGaussPulse(1, 'use', use);
        assert(isa(pulse, 'struct'), ['Expected output to be a struct for use: ', use]);
    end
end
