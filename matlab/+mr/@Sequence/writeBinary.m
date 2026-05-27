function writeBinary(obj,filename,create_signature)
%WRITEBINARY Write sequence to file in binary format.
%   WRITEBINARY(seqObj, filename) Write the sequence data to the given
%   filename using the binary version of the Pulseq open file format for MR
%   sequences. The file specification is available at
%   http://pulseq.github.io
%
%   Examples:
%   Write the sequence file to the sequences directory
%
%       writeBinary(seqObj,'sequences/gre.bseq')
%
% See also  readBinary

if (nargin<3)
    create_signature=true;
end

% handle RequiredExtensions definition (same as write())
if ~isempty(obj.rotationLibrary.keys)
    RD=obj.getDefinition('RequiredExtensions');
    if isempty(RD) || isempty(strfind(RD,'ROTATIONS'))
        RD=mr.aux.strstrip([mr.aux.strstrip(RD) ' ROTATIONS']);
        obj.setDefinition('RequiredExtensions', RD);
    end
end

binaryCodes = obj.getBinaryCodes();
fid=fopen(filename, 'w');
fwrite(fid, binaryCodes.fileHeader, 'int64');
fwrite(fid, int64(obj.version_major), 'int64');
fwrite(fid, int64(obj.version_minor), 'int64');
fwrite(fid, int64(obj.version_revision), 'int64');

if ~isempty(obj.definitions)
    fwrite(fid, binaryCodes.section.definitions, 'int64');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    fwrite(fid, length(keys), 'int64');
    for i = 1:length(keys)
        fwrite(fid, length(keys{i}),'int32');
        fwrite(fid, keys{i},'char');
        val = values{i};
        fwrite(fid, length(val), 'int32');
        if ischar(val)
            fwrite(fid, 'c', 'char');
            fwrite(fid, val, 'char');
        elseif isinteger(val)
            fwrite(fid, 'i', 'char');
            fwrite(fid, val, 'int32');
        elseif isfloat(val)
            fwrite(fid, 'f', 'char');
            fwrite(fid, val, 'float64');
        else
            error(['unknown type of the value type for ' keys{i} ]);
        end
    end
end

% Blocks: write count, then per block: duration (int64, in blockDurationRaster units)
% followed by 6 event IDs (int32): rf, gx, gy, gz, adc, ext
fwrite(fid, binaryCodes.section.blocks, 'int64');
fwrite(fid, length(obj.blockEvents), 'int64');
for i = 1:length(obj.blockEvents)
    bd = obj.blockDurations(i) / obj.blockDurationRaster;
    bdr = round(bd);
    assert(abs(bdr - bd) < 1e-6);
    fwrite(fid, bdr, 'int64');                    % block duration in raster units
    fwrite(fid, obj.blockEvents{i}(2:end), 'int32'); % rf, gx, gy, gz, adc, ext
end

% RF: amp(f64) mag_id(i32) phase_id(i32) time_shape_id(i32) center(i64,us)
%     delay(i64,us) freqPPM(f64) phasePPM(f64) freq(f64) phase(f64) use(char)
% array layout: [amp mag_id phase_id time_shape_id center delay freqPPM phasePPM freq phase]
if ~isempty(obj.rfLibrary.keys)
    keys = obj.rfLibrary.keys;
    fwrite(fid, binaryCodes.section.rf, 'int64');
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.rfLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1), 'float64');              % amp
        fwrite(fid, data(2:4), 'int32');              % mag_id, phase_id, time_shape_id
        fwrite(fid, round(data(5)*1e12), 'int64');    % center (ps)
        fwrite(fid, round(data(6)*1e12), 'int64');    % delay (ps)
        fwrite(fid, data(7:10), 'float64');           % freqPPM, phasePPM, freq, phase
        fwrite(fid, obj.rfLibrary.type(k), 'char');   % use
    end
end

arbGradMask = obj.gradLibrary.type=='g';
trapGradMask = obj.gradLibrary.type=='t';

% Arbitrary gradients: amp(f64) first(f64) last(f64) amp_shape_id(i32)
%                      time_shape_id(i32) delay(i32,us)
% array layout: [amp first last amp_shape_id time_shape_id delay]
if any(arbGradMask)
    keys = obj.gradLibrary.keys;
    fwrite(fid, binaryCodes.section.gradients, 'int64');
    fwrite(fid, length(keys(arbGradMask)), 'int64');
    for k = keys(arbGradMask)
        data = obj.gradLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1:3), 'float64');            % amp, first, last
        fwrite(fid, data(4:5), 'int32');              % amp_shape_id, time_shape_id
        fwrite(fid, round(data(6)*1e12), 'int64');    % delay (ps)
    end
end

% Trapezoid gradients: amp(f64) rise(i64,us) flat(i64,us) fall(i64,us) delay(i64,us)
% array layout: [amp rise flat fall delay]
if any(trapGradMask)
    keys = obj.gradLibrary.keys;
    fwrite(fid, binaryCodes.section.trapezoids, 'int64');
    fwrite(fid, length(keys(trapGradMask)), 'int64');
    for k = keys(trapGradMask)
        data = obj.gradLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1), 'float64');              % amp
        fwrite(fid, data(2:5)*1e12, 'int64');         % rise, flat, fall, delay (ps)
    end
end

% ADC: num(i64) dwell(i64,ns) delay(i64,us) freqPPM(f64) phasePPM(f64)
%      freq(f64) phase(f64) phase_id(i32)
% array layout: [num dwell delay freqPPM phasePPM freq phase phase_id]
if ~isempty(obj.adcLibrary.keys)
    keys = obj.adcLibrary.keys;
    fwrite(fid, binaryCodes.section.adc, 'int64');
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.adcLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1), 'int64');                % num
        fwrite(fid, round(data(2)*1e12), 'int64');    % dwell (ps)
        fwrite(fid, round(data(3)*1e12), 'int64');    % delay (ps)
        fwrite(fid, data(4:7), 'float64');            % freqPPM, phasePPM, freq, phase
        fwrite(fid, data(8), 'int32');                % phase_id
    end
end

if ~isempty(obj.shapeLibrary.keys)
    keys = obj.shapeLibrary.keys;
    fwrite(fid, binaryCodes.section.shapes, 'int64');
    fwrite(fid, length(keys), 'int64');
    for k = keys
        shape = obj.shapeLibrary.data(k).array;
        num_samples = shape(1);
        data = shape(2:end);
        fwrite(fid, k, 'int32');
        fwrite(fid, num_samples, 'int64');            % num uncompressed
        fwrite(fid, length(data), 'int64');           % num compressed
        fwrite(fid, data, 'float32');
    end
end

% Extensions: id(i32) type(i32) ref(i32) next_id(i32)
if ~isempty(obj.extensionLibrary.keys)
    keys = obj.extensionLibrary.keys;
    fwrite(fid, binaryCodes.section.extensions, 'int64');
    fwrite(fid, length(keys), 'int64');
    for k = keys
        fwrite(fid, k, 'int32');
        fwrite(fid, round(obj.extensionLibrary.data(k).array), 'int32'); % type, ref, next_id
    end
end

% Triggers: id(i32) type(i32) channel(i32) delay(i64,ps) duration(i64,ps)
if ~isempty(obj.trigLibrary.keys)
    keys = obj.trigLibrary.keys;
    fwrite(fid, binaryCodes.section.triggers, 'int64');
    fwrite(fid, obj.getExtensionTypeID('TRIGGERS'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.trigLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1:2), 'int32');              % type, channel
        fwrite(fid, round(data(3:4)*1e12), 'int64');   % delay, duration (ps)
    end
end

% Labels (LABELSET and LABELINC): id(i32) value(i32) label_index(i32)
if ~isempty(obj.labelsetLibrary.keys)
    keys = obj.labelsetLibrary.keys;
    fwrite(fid, binaryCodes.section.labelset, 'int64');
    fwrite(fid, obj.getExtensionTypeID('LABELSET'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.labelsetLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1:2), 'int32');              % value, label_index
    end
end

if ~isempty(obj.labelincLibrary.keys)
    keys = obj.labelincLibrary.keys;
    fwrite(fid, binaryCodes.section.labelinc, 'int64');
    fwrite(fid, obj.getExtensionTypeID('LABELINC'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.labelincLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1:2), 'int32');              % value, label_index
    end
end

% Soft delays: id(i32) num(i32) offset(i64,ps) factor(f64) hint_len(i32) hint(chars)
if ~isempty(obj.softDelayLibrary.keys)
    keys = obj.softDelayLibrary.keys;
    fwrite(fid, binaryCodes.section.softdelays, 'int64');
    fwrite(fid, obj.getExtensionTypeID('DELAYS'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        data = obj.softDelayLibrary.data(k).array;
        hint_str = obj.softDelayHints2{data(4)};
        fwrite(fid, k, 'int32');
        fwrite(fid, data(1), 'int32');                % num
        fwrite(fid, round(data(2)*1e12), 'int64');    % offset (ps)
        fwrite(fid, data(3), 'float64');              % factor
        fwrite(fid, length(hint_str), 'int32');       % hint string length
        fwrite(fid, hint_str, 'char');                % hint string (no null terminator)
    end
end

% RF shims: id(i32) num_chan(i32) mag_c1(f64) phase_c1(f64) ...
if ~isempty(obj.rfShimLibrary.keys)
    keys = obj.rfShimLibrary.keys;
    fwrite(fid, binaryCodes.section.rfshims, 'int64');
    fwrite(fid, obj.getExtensionTypeID('RF_SHIMS'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        chan_data = obj.rfShimLibrary.data(k).array;
        fwrite(fid, k, 'int32');
        fwrite(fid, length(chan_data)/2, 'int32');    % num channels
        fwrite(fid, chan_data, 'float64');            % mag/phase pairs
    end
end

% Rotations: id(i32) q0(f64) qx(f64) qy(f64) qz(f64)
if ~isempty(obj.rotationLibrary.keys)
    keys = obj.rotationLibrary.keys;
    fwrite(fid, binaryCodes.section.rotations, 'int64');
    fwrite(fid, obj.getExtensionTypeID('ROTATIONS'), 'int32'); % extension type ID
    fwrite(fid, length(keys), 'int64');
    for k = keys
        fwrite(fid, k, 'int32');
        fwrite(fid, obj.rotationLibrary.data(k).array, 'float64'); % quaternion [q0 qx qy qz]
    end
end

fclose(fid);

if create_signature
    % sign the file (this version with re/loading the file is a factor 2 faster than the sprintf() based one that kept a memory-copy of the data written)

    % re-open and read in the file
    fid=fopen(filename, 'r');
    buf=fread(fid);
    fclose(fid);

    % calculate the digest
    md5hash=mr.aux.md5(buf);
    %fprintf('%s\n',md5hash);

    % store the signature in the seq object
    obj.signatureType='md5';
    obj.signatureFile='bin';
    obj.signatureValue=md5hash;

    % re-open the file for appending
    fid=fopen(filename, 'a');
    fseek(fid, 0, 'eof'); % Octave seems to need this inspite of 'a'
    fpos=ftell(fid);
    fwrite(fid, binaryCodes.section.signature, 'int64');
    % signature type: length,string
    fwrite(fid, length(obj.signatureType), 'int32');
    fwrite(fid, obj.signatureType, 'char');
    % signature: length,data (as bytes, not characters)
    fwrite(fid, length(obj.signatureValue)/2, 'int32');
    for i=1:length(obj.signatureValue)/2
        fwrite(fid, hex2dec(obj.signatureValue(i*2-1:i*2)), 'uint8');
    end
    % the original length of the file prior to adding the signature for easier signature validation : int64
    fwrite(fid,fpos,'int64');
    fclose(fid);
end

end
