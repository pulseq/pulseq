function tests = testCalcDuration
    tests = functiontests(localfunctions);
end

%% Test No Event (Should return duration 0)
function testNoEvent(testCase)
    testCase.verifyEqual(mr.calcDuration([]), 0);
end

%% Define Known Events and Expected Durations
function testSingleEvents(testCase)
    eventZoo = getEventZoo();

    for i = 1:length(eventZoo)
        event = eventZoo(i).event;
        expected_dur = eventZoo(i).expected_dur;
        testCase.verifyEqual(mr.calcDuration(event), expected_dur);
    end
end

%% Test Combinations of 2 Events
function testEventCombinations2(testCase)
    validateEventCombinations(testCase, 2);
end

%% Test Combinations of 3 Events
function testEventCombinations3(testCase)
    validateEventCombinations(testCase, 3);
end

%% Function to generate the known event zoo
function eventZoo = getEventZoo()
    eventZoo(1) = struct('name', 'trapz_amp1', 'event', mr.makeTrapezoid('x', 'amplitude', 1, 'duration', 1), 'expected_dur', 1);
    eventZoo(2) = struct('name', 'trapz_amp1_delayed1', 'event', mr.makeTrapezoid('x', 'amplitude', 1, 'duration', 1, 'delay', 1), 'expected_dur', 2);
    eventZoo(3) = struct('name', 'delay1', 'event', mr.makeDelay(1), 'expected_dur', 1);
    eventZoo(4) = struct('name', 'delay0', 'event', mr.makeDelay(0), 'expected_dur', 0);
    eventZoo(5) = struct('name', 'rf0_block1', 'event', mr.makeBlockPulse(0, 'duration', 1), 'expected_dur', 1);
    eventZoo(6) = struct('name', 'rf10_block1', 'event', mr.makeBlockPulse(10, 'duration', 1), 'expected_dur', 1);
    eventZoo(7) = struct('name', 'rf10_block1_delay1', 'event', mr.makeBlockPulse(10, 'duration', 1, 'delay', 1), 'expected_dur', 2);
    eventZoo(8) = struct('name', 'adc3', 'event', mr.makeAdc(1, 'duration', 3), 'expected_dur', 3);
    eventZoo(9) = struct('name', 'adc3_delayed', 'event', mr.makeAdc(1, 'duration', 3, 'delay', 1), 'expected_dur', 4);
    eventZoo(10) = struct('name', 'outputOsc042', 'event', mr.makeDigitalOutputPulse('osc0', 'duration', 42), 'expected_dur', 42);
    eventZoo(11) = struct('name', 'outputOsc142_delay1', 'event', mr.makeDigitalOutputPulse('osc1', 'duration', 42, 'delay', 1), 'expected_dur', 43);
    eventZoo(12) = struct('name', 'outputExt42_delay9', 'event', mr.makeDigitalOutputPulse('osc1', 'duration', 42, 'delay', 9), 'expected_dur', 51);
    eventZoo(13) = struct('name', 'triggerPhysio159', 'event', mr.makeTrigger('physio1', 'duration', 59), 'expected_dur', 59);
    eventZoo(14) = struct('name', 'triggerPhysio259_delay1', 'event', mr.makeTrigger('physio2', 'duration', 59, 'delay', 1), 'expected_dur', 60);
    eventZoo(15) = struct('name', 'label0', 'event', mr.makeLabel('SET', 'SLC', 0), 'expected_dur', 0);
end


%% Function to validate event combinations
function validateEventCombinations(testCase, numToCombine)
    eventZoo = getEventZoo();
    n = length(eventZoo);
    
    % Generate combinations of numToCombine elements
    combos = nchoosek(1:n, numToCombine); % MATLAB equivalent of itertools.combinations_with_replacement

    for i = 1:size(combos, 1)
        comboIdx = combos(i, :);
        events = {eventZoo(comboIdx).event};
        expected_dur = max([eventZoo(comboIdx).expected_dur]);

        % Verify duration calculation
        testCase.verifyEqual(mr.calcDuration(events{:}), expected_dur);
    end
end
