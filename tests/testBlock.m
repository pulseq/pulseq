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

%% Helper functions to create rotation matrix
function R = rotmat()
    % ROTATION_MATRIX Create a 3x3 rotation matrix for a given angle (degrees)
    theta = deg2rad(90.0);
    R0 = [cos(theta), -sin(theta), 0];
    R1 = [sin(theta),  cos(theta), 0];
    R2 = [0, 0, 1];
    R = mr.makeRotation([R0; R1; R2]);
end

function ID = identity()
    ID = mr.makeRotation(eye(3));
end

%% Tests for gradient continuity in addBlock
function testGradientContinuity1(testCase)
    % Trap followed by extended gradient: No error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap);
    seq.addBlock(testCase.TestData.gx_extended);
    seq.addBlock(testCase.TestData.gx_trap);
end

function testGradientContinuity2(testCase)
    % Trap followed by non-zero gradient: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh));
end

function testGradientContinuity3(testCase)
    % Gradient starts at non-zero in first block: Should raise an error
    seq = mr.Sequence();
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh));
end

function testGradientContinuity4(testCase)
    % Gradient starts and ends at non-zero: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.delay);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_allhigh));
end

function testGradientContinuity5(testCase)
    % Gradient starts at zero and has a delay: No error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_extended_delay);
end

function testGradientContinuity6(testCase)
    % Gradient starts at non-zero in other blocks: Should raise an error
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

%% Tests for gradient continuity in addBlock when blocks are rotated
function testGradientContinuityRot1(testCase)
    % Trap followed by extended gradient: No error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap, identity());
    seq.addBlock(testCase.TestData.gx_extended, identity());
    seq.addBlock(testCase.TestData.gx_trap, identity());
end

function testGradientContinuityRot2(testCase)
    % Trap followed by non-zero gradient: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_trap, identity());
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh, identity()));
end

function testGradientContinuityRot3(testCase)
    % Gradient starts at non-zero in first block: Should raise an error
    seq = mr.Sequence();
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh, identity()));
end

function testGradientContinuityRot4(testCase)
    % Gradient starts and ends at non-zero: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.delay);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_allhigh, identity()));
end

function testGradientContinuityRot5(testCase)
    % Gradient starts at zero and has a delay: No error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_extended_delay, identity());
end

function testGradientContinuityRot6(testCase)
    % Gradient starts at non-zero in other blocks: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.delay);
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh, identity()));
end

function testGradientContinuityRot7(testCase)
    % Non-zero, but non-connecting gradients: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_endshigh, identity());
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh2, identity()));
end

function testGradientContinuityRot8(testCase)
    % Non-zero, both gradients are rotated by the same angle: No error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_endshigh, rotmat());
    seq.addBlock(testCase.TestData.gx_startshigh, rotmat());
end

function testGradientContinuityRot9(testCase)
    % Non-zero, new gradient has different rotation from previous: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_endshigh, identity());
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh, rotmat()));
end

function testGradientContinuityRot10(testCase)
    % Non-zero, new gradient has different rotation from previous: Should raise an error
    seq = mr.Sequence();
    seq.addBlock(testCase.TestData.gx_endshigh, rotmat());
    verifyErrorThrown(testCase, @() seq.addBlock(testCase.TestData.gx_startshigh, identity()));
end


