function write(obj,filename,create_signature)
%WRITE Write sequence to file.
%   WRITE(seqObj, filename) Write the sequence data to the given
%   filename using the open file format for MR sequences.
%
%   Examples:
%   Write the sequence file to the my_sequences directory
%
%       write(seqObj,'my_sequences/gre.seq')
%
% See also  read

if (nargin<3)
    create_signature=true;
end

fid=fopen(filename, 'w');
assert(fid ~= -1, 'Cannot open file: %s', filename);
fprintf(fid, '# Pulseq sequence file\n');
fprintf(fid, '# Created by MATLAB mr toolbox\n\n');

% we always write files in the default current version, which may be
% differen to one, loaded (and stored in the seq object)
[version_major, version_minor, version_revision]=mr.aux.version('output');
fprintf(fid, '[VERSION]\n');
fprintf(fid, 'major %s\n', num2str(version_major));
fprintf(fid, 'minor %s\n', num2str(version_minor));
fprintf(fid, 'revision %s\n', num2str(version_revision));
fprintf(fid, '\n');

% handle RequiredExtensions definition
if ~isempty(obj.rotationLibrary.keys)
    RD=obj.getDefinition('RequiredExtensions');
    if isempty(strfind(RD,'ROTATIONS'))
        RD=mr.aux.strstrip([mr.aux.strstrip(RD) ' ROTATIONS']);
        obj.setDefinition('RequiredExtensions', RD);
    end
end

if ~isempty(obj.definitions)
    fprintf(fid, '[DEFINITIONS]\n');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    for i=1:length(keys)
        fprintf(fid, '%s ', keys{i});
        if (ischar(values{i}))
            fprintf(fid, '%s ', values{i});
        else
            fprintf(fid, '%.9g ', values{i});
        end
        fprintf(fid, '\n');
    end
	fprintf(fid, '\n');
end

fprintf(fid, '# Format of blocks:\n');
fprintf(fid, '# NUM DUR RF  GX  GY  GZ  ADC  EXT\n');
fprintf(fid, '[BLOCKS]\n');
idFormatWidth = length(num2str(length(obj.blockEvents)));
idFormatStr = ['%' num2str(idFormatWidth) 'd'];
for i = 1:length(obj.blockEvents)
    %fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d 0\n'],[i obj.blockEvents(i,:)]);
    %fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d 0\n'],[i obj.blockEvents{i}]);
    bd=obj.blockDurations(i)/obj.blockDurationRaster;
    bdr=round(bd);
    assert(abs(bdr-bd)<1e-6); % this may still trigger false alarms for very long delays due to the limited accuracy of the double 
    fprintf(fid,[idFormatStr ' %3d %3d %3d %3d %3d %2d %2d\n'], ...
            [i bdr obj.blockEvents{i}(2:end)]); 
end
fprintf(fid, '\n');

if ~isempty(obj.rfLibrary.keys)
    fprintf(fid, '# Format of RF events:\n');
    fprintf(fid, '# id ampl. mag_id phase_id time_shape_id center delay freqPPM phasePPM freq phase use\n');
    fprintf(fid, '# ..   Hz      ..       ..            ..     us    us     ppm  rad/MHz   Hz   rad  ..\n');
    fprintf(fid,['# Field ''use'' is the initial of: \n#   ' ...
        strtrim(cell2mat(cellfun(@(x) [x ' '], mr.getSupportedRfUse(), 'UniformOutput', false))) ... 
        '\n']);
    fprintf(fid, '[RF]\n');
    keys = obj.rfLibrary.keys;
    for k = keys
        libData1 = obj.rfLibrary.data(k).array(1:4);
        libData2 = obj.rfLibrary.data(k).array(7:10);
        center = obj.rfLibrary.data(k).array(5)*1e6; % us
        delay = round(obj.rfLibrary.data(k).array(6)/obj.rfRasterTime)*obj.rfRasterTime*1e6; % a bit of a hack: round the delay
        fprintf(fid, '%d %12g %d %d %d %g %g %g %g %g %g %c\n', [k libData1 center delay], libData2, obj.rfLibrary.type(k));
    end
    fprintf(fid, '\n');
end

arbGradMask = obj.gradLibrary.type == 'g';
trapGradMask = obj.gradLibrary.type == 't';

if any(arbGradMask)
    fprintf(fid, '# Format of arbitrary gradients:\n');
    fprintf(fid, '#   time_shape_id of 0 means default timing (stepping with grad_raster starting at 1/2 of grad_raster)\n');    
    fprintf(fid, '# id amplitude first last amp_shape_id time_shape_id delay\n'); 
    fprintf(fid, '# ..      Hz/m  Hz/m Hz/m        ..         ..          us\n');
    fprintf(fid, '[GRADIENTS]\n');
    keys = obj.gradLibrary.keys;
    for k = keys(arbGradMask)
        fprintf(fid, '%d %12g %12g %12g %d %d %d\n', ...
                [k obj.gradLibrary.data(k).array(1:5) ...
                 round(obj.gradLibrary.data(k).array(6)*1e6)]);
    end
    fprintf(fid, '\n');
end

if any(trapGradMask)
    fprintf(fid, '# Format of trapezoid gradients:\n');
    fprintf(fid, '# id amplitude rise flat fall delay\n');
    fprintf(fid, '# ..      Hz/m   us   us   us    us\n');
    fprintf(fid, '[TRAP]\n');
    keys = obj.gradLibrary.keys;
    for k = keys(trapGradMask)
        data = obj.gradLibrary.data(k).array;
        data(2:end) = round(1e6*data(2:end));
        fprintf(fid, '%2d %12g %3d %4d %3d %3d\n', [k data]);
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.adcLibrary.keys)
    fprintf(fid, '# Format of ADC events:\n');
    fprintf(fid, '# id num dwell delay freqPPM phasePPM freq phase phase_id\n');
    fprintf(fid, '# ..  ..    ns    us     ppm  rad/MHz   Hz   rad       ..\n');
    fprintf(fid, '[ADC]\n');
    keys = obj.adcLibrary.keys;
    for k = keys
        data = obj.adcLibrary.data(k).array.*[1 1e9 1e6 1 1 1 1 1];
        fprintf(fid, '%d %d %.0f %.0f %g %g %g %g %d\n', [k data]); 
    end
    fprintf(fid, '\n');
end

%if ~isempty(obj.delayLibrary.keys)
%    fprintf(fid, '# Format of delays:\n');
%    fprintf(fid, '# id delay (us)\n');
%    fprintf(fid, '[DELAYS]\n');
%    keys = obj.delayLibrary.keys;
%    for k = keys
%        fprintf(fid, '%d %d\n', ...
%                [k round(1e6*obj.delayLibrary.data(k).array)]);
%    end
%    fprintf(fid, '\n');
%end

if ~isempty(obj.extensionLibrary.keys)
    fprintf(fid, '# Format of extension lists:\n');
    fprintf(fid, '# id type ref next_id\n');
    fprintf(fid, '# next_id of 0 terminates the list\n');
    fprintf(fid, '# Extension list is followed by extension specifications\n');
    fprintf(fid, '[EXTENSIONS]\n');
    keys = obj.extensionLibrary.keys;
    for k = keys
        fprintf(fid, '%d %d %d %d\n', ...
                [k round(obj.extensionLibrary.data(k).array)]);
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.trigLibrary.keys)
    fprintf(fid, '# Extension specification for digital output and input triggers:\n');
    fprintf(fid, '# id type channel delay (us) duration (us)\n');
%     fprintf(fid, 'extension TRIGGERS 1\n'); % fixme: extension ID 1 is hardcoded here for triggers
    fprintf(fid, ['extension TRIGGERS ',num2str(obj.getExtensionTypeID('TRIGGERS')),'\n']);

    keys = obj.trigLibrary.keys;
    for k = keys
        fprintf(fid, '%d %d %d %d %d\n', ...
                [k round(obj.trigLibrary.data(k).array.*[1 1 1e6 1e6])]); 
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.labelsetLibrary.keys) || ~isempty(obj.labelincLibrary.keys)
    lbls=mr.getSupportedLabels();

    if ~isempty(obj.labelsetLibrary.keys)
       fprintf(fid, '# Extension specification for setting labels:\n');
        fprintf(fid, '# id set labelstring\n');
        tid=obj.getExtensionTypeID('LABELSET');
        fprintf(fid, ['extension LABELSET ',num2str(tid),'\n']);
        keys = obj.labelsetLibrary.keys;
        for k = keys
                fprintf(fid, '%d %d %s\n', ... 
                    k, obj.labelsetLibrary.data(k).array(1),lbls{obj.labelsetLibrary.data(k).array(2)});
        end
        fprintf(fid, '\n');
    end
    if ~isempty(obj.labelincLibrary.keys)
        fprintf(fid, '# Extension specification for increasing labels:\n');
        fprintf(fid, '# id inc labelstring\n');
        tid=obj.getExtensionTypeID('LABELINC');
        fprintf(fid, ['extension LABELINC ',num2str(tid),'\n']);
        lbls=mr.getSupportedLabels();
        keys = obj.labelincLibrary.keys;
        for k = keys
                fprintf(fid, '%d %d %s\n', ... 
                    k, obj.labelincLibrary.data(k).array(1),lbls{obj.labelincLibrary.data(k).array(2)});
        end
        fprintf(fid, '\n');
    end
end

if ~isempty(obj.softDelayLibrary.keys)
    fprintf(fid, '# Extension specification for soft delays:\n');
    fprintf(fid, '# id num offset factor hint\n');
    fprintf(fid, '# ..  ..     us     ..   ..\n');
    fprintf(fid, ['extension DELAYS ',num2str(obj.getExtensionTypeID('DELAYS')),'\n']);

    keys = obj.softDelayLibrary.keys;
    for k = keys
        fprintf(fid, '%d %d %g %g %s\n', ...
                k, obj.softDelayLibrary.data(k).array(1), obj.softDelayLibrary.data(k).array(2)*1e6, obj.softDelayLibrary.data(k).array(3), obj.softDelayHints2{obj.softDelayLibrary.data(k).array(4)});
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.rfShimLibrary.keys)
    fprintf(fid, '# Extension specification for RF shimming:\n');
    fprintf(fid, '# id num_chan magn_c1 phase_c1 magn_c2 phase_c2 ...\n');
    fprintf(fid, ['extension RF_SHIMS ',num2str(obj.getExtensionTypeID('RF_SHIMS')),'\n']);

    keys = obj.rfShimLibrary.keys;
    for k = keys
        fprintf(fid, '%d %d', [k length(obj.rfShimLibrary.data(k).array)/2]); 
        fprintf(fid, ' %g', obj.rfShimLibrary.data(k).array);
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.rotationLibrary.keys)
    fprintf(fid, '# Extension specification for rotation events:\n');
    fprintf(fid, '# id RotQuat0 RotQuatX RotQuatY RotQuatZ\n');
    fprintf(fid, ['extension ROTATIONS ',num2str(obj.getExtensionTypeID('ROTATIONS')),'\n']);

    keys = obj.rotationLibrary.keys;
    for k = keys
        fprintf(fid, '%d ', k );
        fprintf(fid, ' %g', obj.rotationLibrary.data(k).array);
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
end

if ~isempty(obj.shapeLibrary.keys)
    fprintf(fid, '# Sequence Shapes\n');
    fprintf(fid, '[SHAPES]\n\n');
    keys = obj.shapeLibrary.keys;
    for k = keys
        shape_dat = obj.shapeLibrary.data(k).array;
        fprintf(fid, 'shape_id %d\n', k);
        fprintf(fid, 'num_samples %d\n', shape_dat(1));
        fprintf(fid, '%.9g\n', shape_dat(2:end));
        fprintf(fid, '\n');
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
    if mr.aux.isOctave()
      md5hash=hash('MD5',buf); % Octave-specific function
    else
      md5hash=md5_java(buf); % Matlab Java hack
    end
    %fprintf('%s\n',md5hash);
    
    % store the signature in the object
    obj.signatureType='md5';
    obj.signatureFile='text';
    obj.signatureValue=md5hash;
    
    % re-open the file for appending
    fid=fopen(filename, 'a');
    fprintf(fid, '\n[SIGNATURE]\n'); % the preceding new line BELONGS to the signature (and needs to be sripped away to recalculate the signature)
    fprintf(fid, '# This is the hash of the Pulseq file, calculated right before the [SIGNATURE]\n');
    fprintf(fid, '# section was added. It can be reproduced/verified with md5sum if the file\n');
    fprintf(fid, '# trimmed to the position right above [SIGNATURE]. The new line character\n');
    fprintf(fid, '# preceding [SIGNATURE] BELONGS to the signature (and needs to be sripped away\n');
    fprintf(fid, '# for recalculating/verification)\n');
    fprintf(fid, 'Type md5\n');
    fprintf(fid, 'Hash %s\n', md5hash);
    fclose(fid);
end

end

function out=md5_java(buf) 
    import java.security.*;
    import java.math.*;
    import java.lang.String;
    
    md = MessageDigest.getInstance('MD5');
    hash = md.digest(double(buf));
    bi = BigInteger(1, hash);
    
    out=char(String.format('%032x', bi));
end
