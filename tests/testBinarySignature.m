function tests = testBinarySignature
    tests = functiontests(localfunctions);
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

        [sectionCode, sigType, sigHex, signedLen] = readBinarySignatureTrailer(pathBin);
        expectedSectionCode = seqLoaded.getBinaryCodes().section.signature;
        computedHex = computeMd5OfPrefix(pathBin, signedLen);

        testCase.verifyEqual(int64(sectionCode), int64(expectedSectionCode));
        testCase.verifyEqual(lower(sigType), 'md5');
        testCase.verifyEqual(lower(seqLoaded.signatureType), 'md5');
        testCase.verifyEqual(lower(seqLoaded.signatureFile), 'bin');
        testCase.verifyEqual(lower(seqLoaded.signatureValue), lower(sigHex));
        testCase.verifyEqual(lower(sigHex), lower(computedHex));
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
        'timeBwProduct', 4);
    gzReph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1e-3);
    seq2.addBlock(rf, gz);
    seq2.addBlock(gzReph);
    seq2.addBlock(mr.makeDelay(1e-3));

    seqs = {seq1, seq2};
end

function [sectionCode, sigType, sigHex, signedLen] = readBinarySignatureTrailer(filename)
    fid = fopen(filename, 'r');
    assert(fid>=0, 'Failed to open file: %s', filename);
    cleanupObj = onCleanup(@() fclose(fid)); 

    fseek(fid, -8, 'eof');
    signedLen = double(fread(fid, 1, 'int64'));

    fseek(fid, signedLen, 'bof');
    sectionCode = int64(fread(fid, 1, 'int64'));

    typeLen = double(fread(fid, 1, 'int32'));
    sigType = char(fread(fid, typeLen, 'char')');

    hashLen = double(fread(fid, 1, 'int32'));
    hashRaw = uint8(fread(fid, hashLen, 'uint8'));
    if isempty(hashRaw)
        sigHex = '';
    else
        sigHex = lower(reshape(dec2hex(hashRaw, 2)', 1, []));
    end
end

function md5Hex = computeMd5OfPrefix(filename, signedLen)
    fid = fopen(filename, 'r');
    assert(fid>=0, 'Failed to open file: %s', filename);
    cleanupObj = onCleanup(@() fclose(fid)); 

    payload = fread(fid, signedLen, 'uint8=>uint8');
    md5Hex = lower(mr.md5(payload(:)')); % this is a slow but Matlab/Octave compatible function
end

function localCleanup(path)
    if exist(path, 'dir')
        try
            rmdir(path, 's');
        catch
        end
    end
end
