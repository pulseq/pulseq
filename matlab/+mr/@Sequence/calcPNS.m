function [ok, pns_norm, pns_comp, t_axis]=calcPNS(obj,hardware,doPlots)
% calculate PNS using safe model implementation by Szczepankiewicz and Witzel
% assumes safe_pns_prediction package has been downloaded and installed in 
% Matlab path. See http://github.com/filip-szczepankiewicz/safe_pns_prediction
%
% returns pns levels due to respective axes (normalized to 1 and not to 100%)
%
% inputs: 
%     hardware - hardware specifications. see safe_example_hw() from
%                the safe_pns_prediction package. Alternatively a text file
%                in the .asc format (Siemens) can be passed, e.g. for Prisma
%                it is MP_GPA_K2309_2250V_951A_AS82.asc (we leave it as an
%                exercise to the interested user to find were these files
%                can be acquired from);
%     doPlots  - optional parameter (defaluts to true)

if nargin < 3
    doPlots=true;
end

% acquire the entire gradient wave form
gw=obj.waveforms_and_times();
if doPlots
    figure;
    plot(gw{1}(1,:),gw{1}(2,:),gw{2}(1,:),gw{2}(2,:),gw{3}(1,:),gw{3}(2,:)); % plot the entire gradient shape
    title('gradient wave form, in T/m');
end

% find beginning and end times and resample GWs to a regular sampling raster
nt_min=floor(min([gw{1}(1,1) gw{2}(1,1) gw{3}(1,1)])/obj.gradRasterTime+eps); % or should we always start at t=0?
nt_max=ceil(max([gw{1}(1,end) gw{2}(1,end) gw{3}(1,end)])/obj.gradRasterTime-eps);
% shift raster positions to the centers of the raster periods
nt_min = nt_min + 0.5;
nt_max = nt_max - 0.5;
if (nt_min<0.5)
    nt_min=0.5
end
t_axis=(nt_min:nt_max)*obj.gradRasterTime;
gwr=zeros(length(t_axis),3);
for i=1:3
    gwr(:,i)=interp1(gw{i}(1,:),gw{i}(2,:),t_axis,'linear',0);
end

if ischar(hardware)
    % this loads the parameters from the provided text file
    asc=readasc(hardware);
    hardware=asc_to_hw(asc);
end

% use the Szczepankiewicz' and Witzel's implementation
[pns_comp,res]=safe_gwf_to_pns(gwr/obj.sys.gamma, NaN*ones(length(t_axis),1), obj.gradRasterTime, hardware); % the RF vector is unused in the code inside but it is zeropaded and exported ... 
% use the exported RF vector to detect and undo zerpopadding
pns_comp=0.01*pns_comp(~isfinite(res.rf),:)';
% calc pns_norm and the final ok/not_ok
pns_norm=vecnorm(pns_comp);
ok=all(pns_norm<1);
% ready
if doPlots
    % plot results
    figure;safe_plot(pns_comp'*100, obj.gradRasterTime);
end

end

% local utility functions

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


function hw = asc_to_hw(asc)
% function hw = asc_to_hw(asc)
%
% SAFE model parameters for the asc structure as read from the asc file.
% See comments for units.
% 
% Maxim Zaitsev 08/10/2019

if isfield(asc,'asCOMP') && isfield(asc.asCOMP,'tName')
    hw.name          = asc.asCOMP(1).tName;
else
    hw.name          = 'unknown';
end
%hw.look_ahead    =  1.0; % MZ: this is not a real hardware parameter but a coefficient, with which the final result is multiplied

hw.x.tau1        = asc.flGSWDTauX(1);  % ms
hw.x.tau2        = asc.flGSWDTauX(2);  % ms
hw.x.tau3        = asc.flGSWDTauX(3);  % ms
hw.x.a1          = asc.flGSWDAX(1);
hw.x.a2          = asc.flGSWDAX(2);
hw.x.a3          = asc.flGSWDAX(3);
hw.x.stim_limit  = asc.flGSWDStimulationLimitX;  % T/m/s
hw.x.stim_thresh = asc.flGSWDStimulationThresholdX;  % T/m/s
hw.x.g_scale     = asc.asGPAParameters.sGCParameters.flGScaleFactorX;

hw.y.tau1        = asc.flGSWDTauY(1);  % ms
hw.y.tau2        = asc.flGSWDTauY(2);  % ms
hw.y.tau3        = asc.flGSWDTauY(3);  % ms
hw.y.a1          = asc.flGSWDAY(1);
hw.y.a2          = asc.flGSWDAY(2);
hw.y.a3          = asc.flGSWDAY(3);
hw.y.stim_limit  = asc.flGSWDStimulationLimitY;  % T/m/s
hw.y.stim_thresh = asc.flGSWDStimulationThresholdY;  % T/m/s
hw.y.g_scale     = asc.asGPAParameters.sGCParameters.flGScaleFactorY;

hw.z.tau1        = asc.flGSWDTauZ(1);  % ms
hw.z.tau2        = asc.flGSWDTauZ(2);  % ms
hw.z.tau3        = asc.flGSWDTauZ(3);  % ms
hw.z.a1          = asc.flGSWDAZ(1);
hw.z.a2          = asc.flGSWDAZ(2);
hw.z.a3          = asc.flGSWDAZ(3);
hw.z.stim_limit  = asc.flGSWDStimulationLimitZ;  % T/m/s
hw.z.stim_thresh = asc.flGSWDStimulationThresholdZ;  % T/m/s
hw.z.g_scale     = asc.asGPAParameters.sGCParameters.flGScaleFactorZ;
end
