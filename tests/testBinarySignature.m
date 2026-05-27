%!test %%% on Octave run with oruntests() %%%
%! testBinarySignature
function tests = testBinarySignature
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

function test_binary_signature_roundtrip_multiple_sequences(testCase)
    currDir = fileparts(mfilename('fullpath'));

    try
        mr.opts();
    catch
        addpath(fullfile(currDir,'../matlab'));
    end

    seqs = makeTestSequences();

    tmpRoot = tempname;
    mkdir(tmpRoot);
    cleanupObj = onCleanup(@() localCleanup(tmpRoot));

    for iSeq = 1:numel(seqs)
        pathBin = fullfile(tmpRoot, sprintf('signature_%d.bseq', iSeq));

        seq = seqs{iSeq};
        seq.writeBinary(pathBin);

        seqLoaded = mr.Sequence();
        seqLoaded.readBinary(pathBin);

        [signatureValid, storedSignature, computedSignature] = mr.verifyFileSignature(pathBin);

        testCase.verifyEqual(lower(seqLoaded.signatureType), 'md5');
        testCase.verifyEqual(lower(seqLoaded.signatureFile), 'bin');
        testCase.verifyTrue(signatureValid);
        testCase.verifyEqual(storedSignature, computedSignature);
    end
end

function seqs = makeTestSequences()
    seq1 = mr.Sequence();
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    adc = mr.makeAdc(128, 'Duration', 1e-3);
    seq1.addBlock(gx, adc);
    seq1.addBlock(mr.makeDelay(2e-3));

    seq2 = mr.Sequence();
    [rf, gz] = mr.makeSincPulse(pi/2, ...
        'Duration', 2e-3, ...
        'sliceThickness', 5e-3, ...
        'apodization', 0.5, ...
        'timeBwProduct', 4, ...
        'use', 'excitation');
    gzReph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1e-3);
    seq2.addBlock(rf, gz);
    seq2.addBlock(gzReph);
    seq2.addBlock(mr.makeDelay(1e-3));

    seqs = {seq1, seq2};
end

function localCleanup(path)
    if exist(path, 'dir')
        if mr.aux.isOctave()
            confirm_recursive_rmdir(false,'local');
        end
        try
            rmdir(path, 's');
        catch
        end
    end
end
