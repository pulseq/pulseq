% simple test example -- testing ADC objects 
% start it by run(test_mr_makeAdc)
% for more information see see Function-Based unit tests in matlab docu
function tests = testMakeAdc
    tests = functiontests(localfunctions);
end

function test_makeAdc_without_timing_should_fail(testCase)
    fh=@() mr.makeAdc(128);
    verifyError(testCase, fh, ?MException);  
end

function test_makeAdc_given_numSamples_and_Dwell_is_valid(testCase)
    adc=mr.makeAdc(128,'Dwell',10e-6);
    verifyEqual(testCase, adc.numSamples, 128);
    verifyEqual(testCase, adc.dwell, 10e-6);
    verifyTrue(testCase, fieldsValid(adc));
end

function test_makeAdc_given_numSamples_and_Duration_is_valid(testCase)
    adc=mr.makeAdc(128,'Duration',1280e-6);
    verifyEqual(testCase, adc.numSamples, 128);
    verifyEqual(testCase, adc.dwell, 10e-6);
    verifyTrue(testCase, fieldsValid(adc));
end

function ok=fieldsValid(adc)
    fields=sort(fieldnames(adc));
    validFields={'deadTime','delay','dwell','freqOffset','freqPPM','numSamples', 'phaseModulation','phaseOffset','phasePPM','type'};
    ok=length(fields)==length(validFields);
    for i=1:length(fields)
        ok=ok && (strcmp(fields{i},validFields{i}));
    end
    if ~ok
        error(['the following fields are problematic: ' strjoin(setxor(fields,validFields))]);
    end 
end