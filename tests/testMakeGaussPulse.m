%!test %%% on Octave run with oruntests() %%%
%! testMakeGaussPulse
function tests = testMakeGaussPulse
    try
        mr.opts();
    catch
        pulseqPath=fullfile(fileparts(mfilename),'..','matlab');
        addpath(genpath(pulseqPath));
    end
    if exist('functiontests')
        tests = functiontests(localfunctions);
    else
        lf=localfunctions();
        testCase=makeOctaveTestCase();
        for i=1:length(lf)
            f=lf{i};
            n=func2str(f);
            if length(n)>3 && strcmp(n(1:4),'test')
                f(testCase);
                fprintf('Test function %s completed successfully\n', n);
            end
        end
    end
end

function test_make_gauss_pulse(testCase)
    % Test the 'use' parameter for valid and invalid inputs
    try
        mr.makeGaussPulse(1, 'use', 'invalid');
        error('Expected error for invalid use parameter not thrown');
    catch ME
        assert(~isempty(strfind(lower(ME.message), 'of ''use'''))); % Octave and Matlab throw different errors, but both contain "of 'use'"
    end

    % Loop through all supported rf uses and check if make_gauss_pulse works
    supported_rf_uses = mr.getSupportedRfUse(); % Assuming this function returns a list of supported uses

    for i = 1:numel(supported_rf_uses)
        use = supported_rf_uses{i};
        pulse = mr.makeGaussPulse(1, 'use', use);
        assert(isa(pulse, 'struct'), ['Expected output to be a struct for use: ', use]);
    end
end
