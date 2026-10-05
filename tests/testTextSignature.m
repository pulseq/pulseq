%!test %%% on Octave run with oruntests() %%%
%! testTextSignature
function tests = testTextSignature
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

function test_text_signature_roundtrip(testCase)
    seq = mr.Sequence();
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    adc = mr.makeAdc(128, 'Duration', 1e-3);
    seq.addBlock(gx, adc);
    seq.addBlock(mr.makeDelay(2e-3));

    pathSeq = [tempname '.seq'];
    cleanupObj = onCleanup(@() delete(pathSeq));

    seq.write(pathSeq);

    seqLoaded = mr.Sequence();
    seqLoaded.read(pathSeq);

    [signatureValid, storedSignature, computedSignature] = mr.verifyFileSignature(pathSeq);

    testCase.verifyEqual(seq.signatureType, 'md5');
    testCase.verifyEqual(seq.signatureFile, 'text');
    testCase.verifyEqual(seqLoaded.signatureType, seq.signatureType);
    testCase.verifyEqual(seqLoaded.signatureFile, seq.signatureFile);
    testCase.verifyEqual(seqLoaded.signatureValue, seq.signatureValue);
    testCase.verifyTrue(signatureValid);
    testCase.verifyEqual(storedSignature, computedSignature);
end
