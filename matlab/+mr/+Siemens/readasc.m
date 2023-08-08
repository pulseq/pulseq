function [asc, extra] = readasc(filePath);
% reads Siemens ASC ascii-formatted textfile 
% return is a matlab structure with fields from the file.
% usage:
% myAsc= readasc(path,fileName);
% [prot, yaps] = readasc(path,fileName);

% Ralph Strecker 16/02/2005
% Maxim Zaitsev 08/10/2019


%%% read asc file and convert it into a structure
fid= fopen(filePath);
endOfAsc=0;

%nextLine=fgetl(fid); %read next line
%while nextLine~=-1 
while ~feof(fid) 
    openbrack= [];
    closebrack= [];
    nextLine=strtrim(fgetl(fid));
    if strcmp(nextLine,'### ASCCONV END ###') % find end of mrProt in the asc file
        endOfAsc=1;
    end
    if isempty(nextLine) || nextLine(1)=='#' 
        continue;
    end
    indEqualSign= findstr(nextLine,'=');
    if ~isempty(indEqualSign)
        fieldName=deblank(nextLine(1:indEqualSign-1));
        openbrack= findstr(fieldName,'[');
        closebrack= findstr(fieldName,']');
        if ~isempty(openbrack) & ~isempty(closebrack)
            fieldName(openbrack)='(';
            fieldName(closebrack)=')';
%             if strcmp(fieldName(end),fieldName(closebrack))
%                 fieldName(closebrack)='}';
%                 fieldName(openbrack)='{';
%             end
            for k=1:length(openbrack)
                counter= str2num(fieldName(openbrack(k)+1:closebrack(k)-1));
                fieldName= [fieldName(1:openbrack(k)),num2str(counter+1),fieldName(closebrack(k):end)];
            end
        end
        openclosebrack= findstr(fieldName,')(');
        fieldName(openclosebrack)=[];
        fieldName(openclosebrack)=',';
        %if  findstr(fieldName,'atImagedNucleus')
        %    fieldName= [fieldName,'.value'];
        %end

        fieldValue= deblank(nextLine(indEqualSign+2:end));
        com=[strfind(fieldValue,'#') strfind(fieldValue,'//')];
        if ~isempty(com)
            com=min(com);
            fieldValue=fieldValue(1:(com-1));
        end
        
%         if findstr(fieldValue,'i0') || findstr(fieldValue,'i1')
%             ind= findstr(fieldValue,'i');
%             fieldValue=[fieldValue(1:ind),'*',fieldValue(ind+1:end)];
%         end
        if ischar(fieldValue)
            ind=findstr(fieldValue,'"');
            %fieldValue(ind)=[];
            fieldValue(ind)='';
        end            
        
        if length(fieldValue)>1 & strcmp(fieldValue(1:2),'0x') & isempty(findstr(fieldName,'atImagedNucleus'))
            fieldValue= hex2dec(fieldValue(3:end)); %fieldValue is hexadecimal
        elseif ~isempty(str2num(fieldValue))
            fieldValue= str2num(fieldValue);
        end
        
        if ischar(fieldValue) && fieldName(end)==')'
            fieldName(end)='}';
            ib=max(strfind(fieldName,'('));
            fieldName(ib)='{';
        end
        
%         fprintf('< %s > ',nextLine);
%         fprintf('fieldName= %s fieldValue= %s\n', fieldName, mat2str(fieldValue));

        if endOfAsc==0
            eval(['asc.' fieldName '=','fieldValue;']);
        else
            eval(['extra.' fieldName '=','fieldValue;']);
        end
    end
        
end
        
fclose(fid);

end