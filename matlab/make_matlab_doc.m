function make_matlab_doc(mfile,outDir)

try
    if nargin<2
        outDir='html/';
    end
    [filePath,baseName,~] = fileparts(mfile);
    evalCode = isempty(strfind(filePath,'@')) && isempty(strfind(filePath,'+'));
    if evalCode
        addpath(filePath);
    end
    
    publish(mfile,'evalCode',evalCode,'showCode',true,'outputDir',outDir);
    
    %% Fix some things to play nicer with doxygen html formatting
    %continue
    
    htmlFile = [outDir '/' baseName '.html'];
    fid=fopen(htmlFile,'r');
    s=fread(fid,'*char').';
    fclose(fid);
    
    s=regexprep(s,'html,body,div,span','html,body,div.content,span');
    s=regexprep(s,'table th {','table.content th {');
    s=regexprep(s,'table td {','table.content td {');
    %s=regexprep(s,'pre, tt, code {','pre, code {');
    
    fid=fopen(htmlFile,'w');
    fwrite(fid,s);
    
    %% Remove timestamp from images
    f=dir([outDir '/' baseName '*.png']);
    for i=1:length(f)
        I=imread([outDir '/' f(i).name]);
        imwrite(I,[outDir '/' f(i).name],'ImageModTime',1);
    end
catch e
    disp(e)
    fprintf('Error detected. Exiting MATLAB...\n');
    exit();
end
end