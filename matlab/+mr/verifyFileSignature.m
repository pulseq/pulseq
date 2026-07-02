function [res, storedSignature, computedSignature] = verifyFileSignature(pulseqFile)
% verifyFileSignature - verify the signature of the provided Pulseq file.
%   both text and binary Pulseq files are supported. 
%   

    fid = fopen(pulseqFile, 'r');
    assert(fid>=0, 'Failed to open file: %s', pulseqFile);
    cleanupObj = onCleanup(@() fclose(fid)); 
    
    % check whether this is a binary file
    binaryCodes = mr.Sequence.getBinaryCodes();
    magicNum = fread(fid,1,'uint64=>uint64');
    isBinary = magicNum==binaryCodes.fileHeader;

    if isBinary
        [ok, sigType, storedSignature, signedLen] = readBinaryFileSignature(fid,binaryCodes);
    else
        [ok, sigType, storedSignature, signedLen] = readTextFileSignature(fid);
        % need to close and reopen the file to reset the 'text' state of the file handle (stupid Matlab)
        fclose(fid);
        fid = fopen(pulseqFile, 'r');
    end
    
    if ~ok
        error('failed to read the signature');
    end

    computedSignature = computeMd5OfPrefix(fid, signedLen);
    res=strcmp(computedSignature,storedSignature);
end

function [ok, sigType, sigHex, signedLen] = readBinaryFileSignature(fid,binaryCodes)
    fseek(fid, -8, 'eof');
    signedLen = double(fread(fid, 1, 'int64'));

    fseek(fid, signedLen, 'bof');
    sectionCode = int64(fread(fid, 1, 'int64'));
    assert(sectionCode==binaryCodes.section.signature);

    typeLen = double(fread(fid, 1, 'int32'));
    sigType = char(fread(fid, typeLen, 'char')');

    hashLen = double(fread(fid, 1, 'int32'));
    hashRaw = uint8(fread(fid, hashLen, 'uint8'));
    if isempty(hashRaw)
        sigHex = '';
    else
        sigHex = lower(reshape(dec2hex(hashRaw, 2)', 1, []));
    end
    ok = true; % more checks?
end

function [ok, sigType, sigHex, signedLen] = readTextFileSignature(fid)
    sigType='';
    sigHex='';
    signedLen=-1;
    fseek(fid, 0, 'bof');
    fLeadingNewLinePos=-1;
    while ~feof(fid)
        fpos=ftell(fid);
        line = fgetl(fid);
        if signedLen<=0 
            % still looking for the SIGNATURE marker
            if length(line)>=11 && strcmp('[SIGNATURE]',line(1:11))
                assert(fLeadingNewLinePos>0);
                signedLen = fLeadingNewLinePos;
            end
            if length(line)==0
                fLeadingNewLinePos=fpos;
            else
                fLeadingNewLinePos=-1;
            end
        else
            % reading the signature
            if length(line)<=5 || line(1)=='#'
                continue;
            end
            if strcmp('Type ',line(1:5))
                sigType=mr.aux.strstrip(line(6:end));
            end
            if strcmp('Hash ',line(1:5))
                sigHex=lower(mr.aux.strstrip(line(6:end)));
            end
        end
    end
    ok = ~isempty(sigType) && ~isempty(sigHex) && signedLen>0;
end

function md5Hex = computeMd5OfPrefix(fid, signedLen)
    fseek(fid, 0, 'bof');
    payload = fread(fid, signedLen);%, 'uint8=>uint8');
    md5Hex = mr.aux.md5(payload);
end
