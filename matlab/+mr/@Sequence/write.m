function write(obj,filename)
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

fid=fopen(filename,'w');
assert(fid~=-1,'Cannot open file: %s',filename);
fprintf(fid,'# Pulseq sequence file\n');
fprintf(fid,'# Created by MATLAB mr toolbox\n\n');

fprintf(fid, '[VERSION]\n');
fprintf(fid, 'major %s\n', num2str(obj.version_major));
fprintf(fid, 'minor %s\n', num2str(obj.version_minor));
fprintf(fid, 'revision %s\n', num2str(obj.version_revision));
fprintf(fid,'\n');

if ~isempty(obj.definitions)
    fprintf(fid,'[DEFINITIONS]\n');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    for i=1:length(keys)
        fprintf(fid,'%s ',keys{i});
        fprintf(fid,'%g ',values{i});
        fprintf(fid,'\n');
    end
	fprintf(fid,'\n');
end

fprintf(fid,'# Format of blocks:\n');
fprintf(fid,'#  #  D RF  GX  GY  GZ ADC\n');
fprintf(fid,'[BLOCKS]\n');
idFormatWidth=length(num2str(length(obj.blockEvents)));
idFormatStr=['%' num2str(idFormatWidth) 'd'];
for i=1:length(obj.blockEvents)
    %fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d 0\n'],[i obj.blockEvents(i,:)]);
    %fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d 0\n'],[i obj.blockEvents{i}]);
    fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d\n'],[i obj.blockEvents{i}]); % PulSeq standard 1.0 (no control events)
end
fprintf(fid,'\n');

if ~isempty(obj.rfLibrary.keys)
    fprintf(fid,'# Format of RF events:\n');
    fprintf(fid,'# id amplitude mag_id phase_id freq phase\n');
    fprintf(fid,'# ..        Hz   ....     ....   Hz   rad\n');
    fprintf(fid,'[RF]\n');
    keys=obj.rfLibrary.keys;
    for k=keys
        libData=obj.rfLibrary.data(k).array(1:5);
        fprintf(fid,'%d %12g %d %d %g %g\n',[k libData]);
    end
    fprintf(fid,'\n');
end

arbGradMask=obj.gradLibrary.type=='g';
trapGradMask=obj.gradLibrary.type=='t';

if any(arbGradMask)
    fprintf(fid,'# Format of arbitrary gradients:\n');
    fprintf(fid,'# id amplitude shape_id\n');
    fprintf(fid,'# ..      Hz/m     ....\n');
    fprintf(fid,'[GRADIENTS]\n');
    keys=obj.gradLibrary.keys;
    for k=keys(arbGradMask)
        fprintf(fid,'%d %12g %d \n',[k obj.gradLibrary.data(k).array]);
    end
    fprintf(fid,'\n');
end

if any(trapGradMask)
    fprintf(fid,'# Format of trapezoid gradients:\n');
    fprintf(fid,'# id amplitude rise flat fall\n');
    fprintf(fid,'# ..      Hz/m   us   us   us\n');
    fprintf(fid,'[TRAP]\n');
    keys=obj.gradLibrary.keys;
    for k=keys(trapGradMask)
        data=obj.gradLibrary.data(k).array;
        data(2:end)=round(1e6*data(2:end));
        fprintf(fid,'%2d %12g %3d %4d %3d\n',[k data]);
    end
    fprintf(fid,'\n');
end

if ~isempty(obj.adcLibrary.keys)
    fprintf(fid,'# Format of ADC events:\n');
    fprintf(fid,'# id num dwell delay freq phase\n');
    fprintf(fid,'# ..  ..    ns    us   Hz   rad\n');
    fprintf(fid,'[ADC]\n');
    keys=obj.adcLibrary.keys;
    for k=keys
        data=obj.adcLibrary.data(k).array(1:5).*[1 1e9 1e6 1 1];
        fprintf(fid,'%2d %3d %6d %3d %g %g\n',[k data]);
    end
    fprintf(fid,'\n');
end

if ~isempty(obj.delayLibrary.keys)
    fprintf(fid,'# Format of delays:\n');
    fprintf(fid,'# id delay (us)\n');
    fprintf(fid,'[DELAYS]\n');
    keys=obj.delayLibrary.keys;
    for k=keys
        fprintf(fid,'%d %d\n',[k round(1e6*obj.delayLibrary.data(k).array)]);
    end
    fprintf(fid,'\n');
end


if ~isempty(obj.shapeLibrary.keys)
    fprintf(fid,'# Sequence Shapes\n');
    fprintf(fid,'[SHAPES]\n\n');
    keys=obj.shapeLibrary.keys;
    for k=keys
        shape_dat = obj.shapeLibrary.data(k).array;
        fprintf(fid,'shape_id %d\n',k);
        fprintf(fid,'num_samples %d\n',shape_dat(1));
        fprintf(fid,'%g\n',shape_dat(2:end));
        fprintf(fid,'\n');
    end
end

fclose(fid);
end
