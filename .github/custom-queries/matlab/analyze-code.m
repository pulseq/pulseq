% Add the current working directory and subdirectories to the path
addpath(genpath(pwd));

% Find all .m files in the current directory and subdirectories
files = dir('**/*.m');

% Initialize results
results = {};

% Loop through each file and perform basic analysis (customize this part)
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    fid = fopen(filePath, 'r');
    
    if fid == -1
        disp(['Could not open file: ', filePath]);
        continue;
    end
    
    code = fread(fid, '*char')'; % Read file contents
    fclose(fid);
    
    % Custom analysis logic here (example: checking for keywords or patterns)
    issues = {}; % You can add code to detect issues
    
    % Store file path and issues found in results
    if ~isempty(issues)
        results{end+1} = {filePath, issues};
    end
end

% Save results to .mat file
save('code-analysis-results.mat', 'results');
