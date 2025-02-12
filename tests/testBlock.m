function tests = testBlock
    tests = functiontests(localfunctions);
end

%% Setup function to define gradients
function setup(testCase)
    import mr.*

    % Gradient definitions used in tests
    testCase.TestData.gx_trap = makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    testCase.TestData.gx_extended = makeExtendedTrapezoid('x', ...
        'Amplitude', [0, 100000, 0], 'Times', [0, 1e-4, 2e-4]);
    testCase.TestData.gx_extended_delay = makeExtendedTrapezoid('x', ...
        'Amplitude', [0, 100000, 0], 'Times', [1e-4, 2e-4, 3e-4]);
    testCase.TestData.gx_endshigh = makeExtendedTrapezoid('x', ...
        'Amplitude', [0, 100000, 100000], 'Times', [0, 1e-4, 2e-4]);
    testCase.TestData.gx_startshigh = makeExtendedTrapezoid('x', ...
        'Amplitude', [100000, 100000, 0], 'Times', [0, 1e-4, 2e-4]);
    testCase.TestData.gx_startshigh2 = makeExtendedTrapezoid('x', ...
        'Amplitude', [200000, 100000, 0], 'Times', [0, 1e-4, 2e-4]);
    testCase.TestData.gx_allhigh = makeExtendedTrapezoid('x', ...
        'Amplitude', [100000, 100000, 100000], 'Times', [0, 1e-4, 2e-4]);
    testCase.TestData.delay = makeDelay(1e-3);

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

%% Tests for gradient continuity in addBlock
function testGradientContinuity1(testCase)
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap);
    seq.addBlock(testCase.TestData.gx_extended);
    seq.addBlock(testCase.TestData.gx_trap);
end

function testGradientContinuity2(testCase)
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh));
end

function testGradientContinuity3(testCase)
    seq = mr.Sequence();
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh));
end

function testGradientContinuity4(testCase)
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.delay);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_allhigh));
end

function testGradientContinuity5(testCase)
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_extended_delay);
end

function testGradientContinuity6(testCase)
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.delay);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh));
end

function testGradientContinuity7(testCase)
    % Non-zero, but non-connecting gradients
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_endshigh);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh2));
end
