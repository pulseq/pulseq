function tests = testCalcADCSegments
    tests = functiontests(localfunctions);
end

%% Setup function: Load test data
function setup(testCase)
    import mr.*
    dirpath = fileparts(mfilename('fullpath')); % Get script directory

    % Load expected output data
    data = readmatrix(fullfile(dirpath, 'expected_output', 'pulseq_calcAdcSeg.txt'));

    % Assign structured data fields
    testCase.TestData.dwell = data(:, 1);
    testCase.TestData.num_samples = data(:, 2);
    testCase.TestData.adc_limit = data(:, 3);
    testCase.TestData.adc_divisor = data(:, 4);
    testCase.TestData.mode = data(:, 5);
    testCase.TestData.res_num_seg = data(:, 6);
    testCase.TestData.res_num_samples_seg = data(:, 7);
    
    % Define system options
    testCase.TestData.system = opts;
    testCase.TestData.system.adcRasterTime = 1e-7;
    testCase.TestData.system.gradRasterTime = 1e-5;
end

%% Helper function to check if an error is thrown
function verifyErrorThrown(testCase, funcHandle)
    didError = false;
    try
        funcHandle();
    catch
        didError = true;
    end
    testCase.verifyTrue(didError, 'Expected an error, but none was thrown.');
end

%% Main Test Function
function testCalcADC(testCase)
    sys = testCase.TestData.system;
    
    % Loop through all test cases
    for i = 1:length(testCase.TestData.dwell)
        dwell = testCase.TestData.dwell(i);
        num_samples = testCase.TestData.num_samples(i);
        adc_limit = testCase.TestData.adc_limit(i);
        adc_divisor = testCase.TestData.adc_divisor(i);
        res_num_seg = testCase.TestData.res_num_seg(i);
        res_num_samples_seg = testCase.TestData.res_num_samples_seg(i);
        modeNum = testCase.TestData.mode(i);
        
        % Convert mode number to string
        if modeNum == 1
            mode = 'shorten';
        else
            mode = 'lenghen';
        end

        % Update system settings
        sys.adcSamplesLimit = adc_limit;
        sys.adcSamplesDivisor = adc_divisor;

        % Run function under test
        [num_seg, num_samples_seg] = mr.calcAdcSeg(num_samples, dwell, sys, mode);
        
        % Validate results against expected values
        testCase.verifyEqual(num_seg, res_num_seg, 'AbsTol', 1e-9);
        testCase.verifyEqual(num_samples_seg, res_num_samples_seg, 'AbsTol', 1e-9);

        % Check if segment samples are within the sample limit
        testCase.verifyLessThanOrEqual(num_samples_seg, adc_limit);

        % Compute durations
        seg_duration = num_samples_seg * dwell;
        adc_duration = seg_duration * num_seg;

        % Check segment alignment with grad raster time
        testCase.verifyTrue(abs(round(seg_duration / sys.gradRasterTime) - (seg_duration / sys.gradRasterTime)) < 1e-9, ...
            'Segment duration is not aligned with grad raster time');

        % Check total ADC alignment with grad raster time
        testCase.verifyTrue(abs(round(adc_duration / sys.gradRasterTime) - (adc_duration / sys.gradRasterTime)) < 1e-9, ...
            'ADC duration is not aligned with grad raster time');

        % Display progress
        if mod(i, 1000) == 0 || i == length(testCase.TestData.dwell)
            fprintf('Progress: %.2f%%\n', i / length(testCase.TestData.dwell) * 100);
        end
    end
end
