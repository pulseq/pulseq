function testCase = makeOctaveTestCase
% Minimal testCase-compatible verifier for Octave fallback execution.

testCase.verifyTrue = @(cond, varargin) iAssertTrue(cond, varargin{:});
testCase.verifyFalse = @(cond, varargin) iAssertFalse(cond, varargin{:});
testCase.verifyEqual = @(a, b, varargin) iAssertEqual(a, b, varargin{:});
testCase.verifyGreaterThan = @(a, b, varargin) iAssertGreaterThan(a, b, varargin{:});
testCase.verifyGreaterThanOrEqual = @(a, b, varargin) iAssertGreaterThanOrEqual(a, b, varargin{:});
testCase.verifyLessThan = @(a, b, varargin) iAssertLessThan(a, b, varargin{:});
testCase.verifyLessThanOrEqual = @(a, b, varargin) iAssertLessThanOrEqual(a, b, varargin{:});
testCase.verifyNotEmpty = @(a, varargin) iAssertNotEmpty(a, varargin{:});
testCase.assertNotEmpty = @(a, varargin) iAssertNotEmpty(a, varargin{:});
testCase.verifyFail = @(varargin) iAssertFail(varargin{:});
testCase.verifyError = @(funcHandle, varargin) iAssertError(funcHandle, varargin{:});
testCase.verifyWarning = @(funcHandle, varargin) iAssertWarning(funcHandle, varargin{:});
testCase.assumeTrue = @(cond, varargin) iAssumeTrue(cond, varargin{:});
testCase.log = @(varargin) fprintf(varargin{:});

end

function iAssertTrue(cond, varargin)
assert(all(cond(:)), iMsg(varargin, 'Expected condition to be true.'));
end

function iAssumeTrue(cond, varargin)
  if ~all(cond(:))
    iMsg(varargin, 'Expected condition to be true, the following test cases are expected to fail.');
  end
end

function iAssertFalse(cond, varargin)
assert(~any(cond(:)), iMsg(varargin, 'Expected condition to be false.'));
end

function iAssertEqual(a, b, varargin)
    if length(varargin)>1 && strcmp('AbsTol',varargin{1})
        assert(abs(a-b)<=varargin{2}, iMsg(varargin(3:end), 'Values are not equal up to the given tolerance.'));
    elseif length(varargin)>1 && strcmp('RelTol',varargin{1})
        assert(abs(a-b)/(0.5*abs(a)+0.5*abs(b))<=varargin{2}, iMsg(varargin(3:end), 'Values are not equal up to the given relative tolerance.'));
    else
        assert(isequaln(a, b), iMsg(varargin, 'Values are not equal.'));
    end
end

function iAssertGreaterThan(a, b, varargin)
assert(all(a(:) > b(:)), iMsg(varargin, 'Expected first value to be greater than second.'));
end

function iAssertGreaterThanOrEqual(a, b, varargin)
assert(all(a(:) >= b(:)), iMsg(varargin, 'Expected first value to be greater than or equal to second.'));
end

function iAssertLessThan(a, b, varargin)
assert(all(a(:) < b(:)), iMsg(varargin, 'Expected first value to be less than second.'));
end

function iAssertLessThanOrEqual(a, b, varargin)
assert(all(a(:) <= b(:)), iMsg(varargin, 'Expected first value to be less than or equal to second.'));
end

function iAssertNotEmpty(a, varargin)
assert(~isempty(a), iMsg(varargin, 'Expected value to be non-empty.'));
end

function iAssertFail(varargin)
assert(false, iMsg(varargin, 'Test failed.'));
end

function iAssertError(funcHandle, varargin)
expectedId = '';
if ~isempty(varargin) && ischar(varargin{1})
    expectedId = varargin{1};
end
caught = false;
caughtId = '';
try
    funcHandle();
catch err
    caught = true;
    if isstruct(err) && isfield(err, 'identifier')
        caughtId = err.identifier;
    end
end
assert(caught, iMsg(varargin, 'Expected an error, but no error was thrown.'));
if ~isempty(expectedId)
    assert(strcmp(caughtId, expectedId), sprintf('Expected error id %s but got %s.', expectedId, caughtId));
end
end

function iAssertWarning(funcHandle, varargin)
expectedId = '';
if ~isempty(varargin) && ischar(varargin{1})
    expectedId = varargin{1};
end
lastwarn('');
funcHandle();
[~, warnId] = lastwarn();
assert(~isempty(warnId), iMsg(varargin, 'Expected warning, but no warning was raised.'));
if ~isempty(expectedId)
    assert(strcmp(warnId, expectedId), sprintf('Expected warning id %s but got %s.', expectedId, warnId));
end
end

function msg = iMsg(args, defaultMsg)
msg = defaultMsg;
for i = numel(args):-1:1
    if ischar(args{i})
        msg = args{i};
        return;
    end
end
end
