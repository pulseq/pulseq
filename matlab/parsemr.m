function parsemr(seqFiles)
%PARSEMR Load a sequence file and display a summary of the sequence events
%   PARSEMR(filename) Load a single file and display a summary of the
%   sequence
%
%   PARSEMR({file1,file2}) Load each filename in the cell array and display
%   a summary of the sequence
%
%   See also mr.Sequence

try
    seq=mr.Sequence();
    seqFiles = cellstr(seqFiles);
    
    for i=1:length(seqFiles)
        file=seqFiles{i};
        
        try
            % Read text file
            seq.read(file,'detectRFuse');
        catch
            % Read binary file
            seq.readBinary(file);
        end
        
        % actually analyze the sequence
        report=seq.testReport;
        
        % Open output file for writing
        [outDir,outFile,~]=fileparts(file);
        if ~isempty(outDir) 
            outDir=[outDir '/'];
        end
        fid=fopen([outDir outFile '.matlab.out'],'w');
        fprintf(fid,'Testing file: %s\n',file);
        fprintf(fid,'%s', [report{:}]);
        fclose(fid);
    end
catch e
    disp(e)
    fprintf('Error detected, exiting MATLAB...\n');
    rethrow(e)
end
