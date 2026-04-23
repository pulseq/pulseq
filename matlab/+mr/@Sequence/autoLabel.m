function [labels, aux] = autoLabel(seq, varargin)
% Automatically find values of the entire sequence or its part.
%
%   autoLabel(seqObj) automatically detects the label evolution from the
%   k-space analysis and applies it to the sequence. Return value of the
%   function is the structure 'labels' with fields named after the labels
%   used in the sequence containing the evolution of the labes for alls
%   ADCs. This structure is similar to that returned by evalLabels().
%
%   autoLabel(...,'blockRange',[first last]) Evaluate label
%   values starting from the first specified block to the last
%   one.
%
%   autoLabel(...,'reflect',[1,2]) Reflect k-space trajectories along
%   directions 1 and 2 (any number of directions between 0 and 3 can be
%   specified) prior to generating label ranges. If used in combination
%   with 'reorder', it is performed first.
%
%   autoLabel(...,'reorder',[2,1]) Reorder axes of k-space trajectories
%   (in this example by swapping directions 1 and 2 prior to generating
%   label ranges (reordering vector can contain 2 or 3 entries). If used in
%   combination with 'reflsect', it is performed last.
%
%   autoLabel(...,'useLabels',labels_struct) Skip label evolution
%   detection and only apply externally-calculated labels to the sequence.
%
%   autoLabel(...,'useAux',aux_struct) use externally-calculated aux struct
%   to apply definitions the sequence. Une combination with 'useLabels'.
%
%   autoLabel(...,'skipApply',true) Skip applying label evolution to the
%   sequence and onoly return the labaels struct.
%
%   see evalLabels()

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'autoLabel';
    parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
    parser.addParamValue('useLabels',struct([]),@(x)(isempty(x) || isstruct(x)));
    parser.addParamValue('useAux',struct([]),@(x)(isempty(x) || isstruct(x)));    
    parser.addParamValue('skipApply',false,@(x)(islogical(x) && isscalar(x)));
    parser.addParamValue('reflect', [], @(x)(isnumeric(x) && numel(x)<4));
    parser.addParamValue('reorder', [], @(x)(isnumeric(x) && numel(x)<4));
end
parse(parser,varargin{:});
opt = parser.Results;

if ~isempty(opt.useLabels) && (~isempty(opt.reflect) || ~isempty(opt.reorder))
    error('Optional parameters ''reflect'' or ''reorder'' only effective for the detection part and cannot be used together with ''useLabels''');
end

if ~isempty(opt.reflect) && numel(opt.reflect)~=numel(unique(opt.reflect))
    error('All indices in ''reflect'' must be unique');
end

if ~isempty(opt.reorder)
    if numel(opt.reorder)~=numel(unique(opt.reorder))
        error('All indices in ''reorder'' must be unique');
    end
    if numel(opt.reorder)~=3 && min(opt.reorder)==1 && max(opt.reorder)==2
        error('If ''reorder'' contains two indices they must be [1 2] or [2 1]');
    end
end

if ~isfinite(opt.blockRange(2))
    opt.blockRange(2)=length(seq.blockEvents);
end

aux = struct();

%% index ADCs, part 1
blockStartTimes=cumsum([0 seq.blockDurations]);
blockStartTimes=blockStartTimes(opt.blockRange(1):opt.blockRange(2));
adcLengths=[];
b_adc=[];
for iB=opt.blockRange(1):opt.blockRange(2)
    block = seq.getBlock(iB);
    if ~isempty(block.adc)
        adcLengths(end+1)=block.adc.numSamples;
        b_adc(end+1)=iB-opt.blockRange(1)+1;
    end
end
adcStartCounts=cumsum([1, adcLengths(1:end-1)]);
% b_adc=zeros(size(t_adcStarts));
% for i=1:numel(t_adcStarts)
%     b_adc(i)=find(blockStartTimes<t_adcStarts(i),1,'last');
% end

if isempty(opt.useLabels)

    %% preparing calculations
    [ktraj_adc, t_adc, ~, ~, t_excitation, ~, slicepos, t_slicepos, gw_pp] = seq.calculateKspacePP('blockRange',opt.blockRange);
    
    t_adcStarts=t_adc(adcStartCounts);
    firstNonNoiseAdc=find(t_adcStarts>t_excitation(1),1,'first');
    firstNonNoiseAdcSampe=find(t_adc>t_excitation(1),1,'first');
    nReadouts = numel(t_adcStarts)-firstNonNoiseAdc+1;
    
    sliceGrads=zeros(size(slicepos));
    for i=1:3
        sliceGrads(i,:)=ppval(gw_pp{i},t_slicepos);
    end
    
    if ~isempty(opt.reflect)
        ktraj_adc(opt.reflect,:)=-ktraj_adc(opt.reflect,:);
        slicepos(opt.reflect,:)=-slicepos(opt.reflect,:);
        sliceGrads(opt.reflect,:)=-sliceGrads(opt.reflect,:);
    end
    if ~isempty(opt.reorder)
        ktraj_adc(1:numel(opt.reorder),:)=ktraj_adc(opt.reorder,:);
        slicepos(1:numel(opt.reorder),:)=slicepos(opt.reorder,:);
        sliceGrads(1:numel(opt.reorder),:)=sliceGrads(opt.reorder,:);
    end
    %% slice positions
    [~,i]=max(abs(sliceGrads));
    mainSliceGradSigns=sign(sliceGrads(sub2ind(size(sliceGrads),i,1:size(sliceGrads,2))));
    sliceNormals=normalize(sliceGrads,'norm').*mainSliceGradSigns(ones(1,3),:);
    sliceOffsets=dot(slicepos,sliceNormals);
    sliceOffsets(~isfinite(sliceOffsets))=0;
    [uniqueSlicePositions, I, sliceCountersAcquisitionOrder] = unique(sliceOffsets,'stable');
    t_usp = t_slicepos(I);
    b_usp = zeros(size(t_usp));
    for i=1:numel(t_usp)
        b_usp(i) = find(blockStartTimes<t_usp(i),1,'last');
    end
    [sortedSlicePositions, ~, sliceCountersSorted] = unique(sliceOffsets);    

    % aux.SliceThickness
    % aux.SliceGap
    % todo: more checks, e.g. if all slice thickeness and gaps are the same...
    GAs1=vecnorm(sliceGrads(:,1));
    if GAs1>0
        bSlice1=find(blockStartTimes<t_slicepos(1),1,'last');
        b=seq.getBlock(bSlice1);
        bw=mr.calcRfBandwidth(b.rf);
        aux.sliceThickness=bw/GAs1;
        if length(sortedSlicePositions)>1
            aux.sliceGap=diff(sortedSlicePositions(1:2))-aux.sliceThickness;
        end
    end

    % useful results:
    %   uniqueSlicePositions : speaks for itself, sorted in the order of occurence
    %   sliceCountersAcquisitionOrder : slice cointers in the acquisition order
    %   b_usp: block indices containg RF pulses corresponding to the first
    %          occurence of each of the unique slice positions
    %   sortedSlicePositions : slice positions in the accending order
    %   sliceCountersSorted : slice counters corresponding to the sorted slice order
    
    %% index ADCs, part 2
    sliceCountersAdc=zeros(1,nReadouts);
    for i=firstNonNoiseAdc:numel(t_adcStarts)
        j=find(t_slicepos<t_adcStarts(i),1,'last');
        sliceCountersAdc(i-firstNonNoiseAdc+1)=find(uniqueSlicePositions==sliceOffsets(j),1,'first');
    end
    % useful results: 
    %   b_adc : vector containing block indices with ADC objects
    %   adcLengths : numbers of samples of each of the ADC objects
    %   t_adcStarts : time point of each ADCs first sample
    %   firstNonNoiseAdc : index of the first ADC after the noise scan (1 means no noise scan)
    %   firstNonNoiseAdcSample : index of the first ADC sample after the noise scan
    %   sliceCountersAdc : slice counters per ADC
    
    %% readout analysis    
    cEchoPos = zeros(1,nReadouts);
    gradReadout = zeros(3,nReadouts);
    signReadout = zeros(1,nReadouts);
    [~,kspaceCenterSample]=min(vecnorm(ktraj_adc(:,firstNonNoiseAdcSampe:end)));
    kspaceCenterPoint=ktraj_adc(:,kspaceCenterSample);
    isCartesianReadout=true;
    for i=1:nReadouts
        c1=adcStartCounts(i+firstNonNoiseAdc-1);
        c2=c1+adcLengths(i+firstNonNoiseAdc-1)-1;
        % find echo pos for each readout
        [~,cEchoPos(i)]=min(vecnorm(ktraj_adc(:,c1:c2)-kspaceCenterPoint)); % for echos positioned between samples there can be some jitter... to reduce jitter we compare not to 0 but to the smallest absolute...
        for j=1:3
            gradReadout(j,i) = ppval(gw_pp{j},t_adc(c1+cEchoPos(i)-1));        
        end
        if ~isempty(opt.reflect)
            gradReadout(opt.reflect,i)=-gradReadout(opt.reflect,i);
        end
        if ~isempty(opt.reorder)
            gradReadout(1:numel(opt.reorder),i)=gradReadout(opt.reorder,i);
        end
        % store the most central readout trajectory as a reference
        if c1<=kspaceCenterSample && c2>=kspaceCenterSample
            cCentralReadout=i; % counter corresponding to the central readout
            kCentralReadout=ktraj_adc(:,c1:c2);
            [~,cCentralReadoutCenter]=min(vecnorm(kCentralReadout));
        end
        % detect readout alternation like for EPI
        [~,j]=max(abs(gradReadout(:,i)));
        signReadout(i)=sign(gradReadout(j,i));
        % detect if the readout is sufficiently the same to see if the trajectory is Cartesian-like
        if i>1 && vecnorm(gradReadout(:,i)*signReadout(i)-gradReadout(:,1)*signReadout(1))>1e-6 % why such threshold? 
            isCartesianReadout=false;
        end
        % % dictionary of readouts
        % kReadoutKey = {ktraj_adc(:,c1:c2) - ktraj_adc(:,c1+cEchoPos(i)-1)};
        % if i==1
        %     readoutDict=dictionary(kReadout,i);
        % else
        %     if ~readoutDict.isKey(kReadout)
        %         readoutDict(kReadout)=i;
        %     end
        % end
    end

    % more readout analysis    
    kCentralReadoutDirection=normalize(gradReadout(:,cCentralReadout),'norm');
    kCentralReadoutProjection=dot(kCentralReadout,kCentralReadoutDirection(:,ones(1,size(kCentralReadout,2))));
    dkR=median(diff(kCentralReadoutProjection));

    % bipolar gradients: detect both and find the reflection point
    if any(diff(signReadout)~=0)
        % find an inverted sign readout close to k-space center
        adcCenters=adcStartCounts(firstNonNoiseAdc:end)+cEchoPos-1;
        isInverted=signReadout==-signReadout(cCentralReadout);
        [~,cInvertedReadout]=min(vecnorm(ktraj_adc(:,adcCenters(isInverted))));
        % here cInvertedReadout is the counter only amongst the inverted readouts, we need to restore its meaning for the total readouts vector
        countInverted=cumsum(isInverted);
        cInvertedReadout=find(cInvertedReadout==countInverted,1);
        assert(isInverted(cInvertedReadout));
        %
        c1=adcStartCounts(firstNonNoiseAdc+cInvertedReadout-1);
        c2=c1+adcLengths(firstNonNoiseAdc+firstNonNoiseAdc-1)-1;
        kInvertedReadoutProjection=dot(ktraj_adc(:,c1:c2),kCentralReadoutDirection(:,ones(1,size(kCentralReadout,2))));
        if signReadout(cCentralReadout)==1
            aux.kReadout = [kCentralReadoutProjection; kInvertedReadoutProjection]; % 1: positive readout, 2: negative readout
        else
            aux.kReadout = [kInvertedReadoutProjection; kCentralReadoutProjection]; % 1: positive readout, 2: negative readout
        end
        figure; plot(aux.kReadout(1,:),'.-'); hold on; plot(aux.kReadout(2,:),'.-');plot(aux.kReadout(2,end:-1:1),'x'); title('bipolar readout alignment'); % if crosses overlap the line then REV is symmetric
    else
        aux.kReadout = kCentralReadoutProjection;
    end
    % not needed information % aux.signReadout=signReadout(cCentralReadout);

    % useful results: 
    %   isCartesianReadout : is this trajectory sufficiently Cartesian-like
    %   gradReadout : readout gradients as 3D vectors
    %   signReadout : sign of the dominant readout component
    %   kCentralReadout : the most central readout
    %   cCentralReadoutCenter : center of the above
    %   cEchoPos : position of the echo for each readout
    
    %% Accelerate calculation for Cartesian-like sequences by only considering the echo position
    if isCartesianReadout
        ktraj_adc = ktraj_adc(:, adcStartCounts(firstNonNoiseAdc:end)+cEchoPos-1);
    else
        error('autoLabel() only supports sufficiently Cartesian sequences');
    end
    
    %% Analyze the trajectory data 
    
    dkR=median(diff(kCentralReadoutProjection));
    
    k_extent=max(abs(ktraj_adc-kspaceCenterPoint),[],2);
    k_scale=max(k_extent);
    k_threshold=abs(dkR/50);
    
    % detect unused dimensions and delete them
    if any(k_extent<k_threshold)
        ktraj_adc(k_extent<k_threshold,:)=[]; % delete rows
        kspaceCenterPoint(k_extent<k_threshold)=[];
        k_extent(k_extent<k_threshold)=[];
    end
    
    % detect dK, k-space reordering and repetitions (or slices, etc)
    kt_sorted=sort(ktraj_adc-kspaceCenterPoint,2);
    dk_all=kt_sorted(:,2:end)-kt_sorted(:,1:(end-1));
    dk_all(dk_all<k_threshold)=NaN;
    dk_min=min(dk_all,[],2);
    dk_max=max(dk_all,[],2);
    dk_all(dk_all-dk_min(:,ones(1,size(dk_all,2)))>k_threshold)=NaN;
    dk_all_cnt=sum(isfinite(dk_all),2);
    dk_all(~isfinite(dk_all))=0;
    dk=sum(dk_all,2)./dk_all_cnt;
    dk(~isfinite(dk))=0; % dk vector ready
    [~,k0_ind]=min(sum(ktraj_adc.^2,1)); % k-space center
    kindex=round((ktraj_adc-ktraj_adc(:,k0_ind*ones(1,size(ktraj_adc,2))))./dk(:,ones(1,size(ktraj_adc,2))));
    kindex(~isfinite(kindex))=0;
    kindex_min=min(kindex,[],2);
    kindex_mat=kindex-kindex_min(:,ones(1,size(ktraj_adc,2)))+1;
    kindex_end=max(kindex_mat,[],2);
    sampler=zeros(length(uniqueSlicePositions),prod(kindex_end));
    repeat=zeros(1,size(ktraj_adc,2));
    for i=1:size(kindex_mat,2)
        switch size(kindex_mat,1)
            case 3
                ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i),kindex_mat(3,i));
            case 2
                ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i)); 
            otherwise
                ind=kindex_mat(1,i); %sub2ind(kindex_end,kindex_mat(1,i));
        end
        r=sampler(sliceCountersAdc(i),ind);
        repeat(i)=r;
        sampler(sliceCountersAdc(i),ind)=r+1;
    end
    % if (max(repeat(:))>0)
    %     kindex=[kindex;(repeat+1)];
    %     kindex_mat=[kindex_mat;(repeat+1)];
    %     kindex_end=max(kindex_mat,[],2);
    % end
    %figure; plot(kindex(1,:),kindex(2,:),'.-');
    %figure; plot(kindex_mat(1,:),'.-');

    %% create labels struct
    labels=struct();
    % NOISE
    % SLC
    % REV
    % LIN,
    % PAR
    % REP (in progress)
    
    % TODO: SEG,ECO,REF,IMA,AVG,SET

    nADCs=nReadouts+firstNonNoiseAdc-1;
    if firstNonNoiseAdc>1
        labels.NOISE=zeros(1,nADCs);
        labels.NOISE(1:(firstNonNoiseAdc-1))=1;
    end
    if any(sliceCountersAdc~=1)
        if firstNonNoiseAdc>1
            labels.SLC=zeros(1,nADCs);
            labels.SLC(firstNonNoiseAdc:end)=sliceCountersAdc-1;
        else
            labels.SLC=sliceCountersAdc-1;
        end
    end
    if any(signReadout<0)
        if firstNonNoiseAdc>1
            labels.REV=zeros(1,nADCs);
            labels.REV(firstNonNoiseAdc:end)=(signReadout<0);
        else
            labels.REV=(signReadout<0);
        end
    end
    if any(kindex_mat(1,:)~=1)
        if firstNonNoiseAdc>1
            labels.LIN=zeros(1,nADCs);
            labels.LIN(firstNonNoiseAdc:end)=kindex_mat(1,:)-1;
        else
            labels.LIN=kindex_mat(1,:)-1;
        end
    end
    if size(kindex_mat,1)>1 && any(kindex_mat(2,:)~=1)
        if firstNonNoiseAdc>1
            labels.PAR=zeros(1,nADCs);
            labels.PAR(firstNonNoiseAdc:end)=kindex_mat(2,:)-1;
        else
            labels.PAR=kindex_mat(2,:)-1;
        end
    end
    if max(repeat(:))>0 
        if firstNonNoiseAdc>1
            labels.REP=zeros(1,nADCs);
            labels.REP(firstNonNoiseAdc:end)=repeat(:).';
        else
            labels.REP=repeat(:).';
        end
    end

    if isfield(labels,'LIN'), aux.kSpaceCenterLine = labels.LIN(cCentralReadout); end
    if isfield(labels,'PAR'), aux.kSpaceCenterPartition = labels.PAR(cCentralReadout); end
    aux.kSpaceCenterSample = cCentralReadoutCenter-1;
    if length(uniqueSlicePositions)>1, aux.SlicePositions = uniqueSlicePositions; end
    % aux.kSpacePhaseEncodingLines
    % aux.PhaseResolution
    % aux.ReadoutOversamplingFactor
    % aux.SliceGap
    % aux.SlicePositions
    % aux.SliceThickness
    % aux.TargetGriddedSamples
    % aux.TrapezoidGriddingParameterss
    % aux.AccelerationFactorPE
    % aux.AccelerationFactor3D
    % aux.FirstFourierLine
    % aux.FirstRefLine
    % aux.FirstFourier3D
    % aux.FirstRef3D

else
    labels=opt.useLabels;
end

if ~isempty(opt.useAux)
    aux=opt.useAux;
end

%% apply labels

if ~opt.skipApply
    lblNames=fieldnames(labels);
    warnBcp = warning ('off','mr:fixmePreviousRotationExtension');
    for i=1:numel(b_adc)
        blkLabels={};
        
        for j=1:length(lblNames)
            if i==1 && labels.(lblNames{j})(1)~=0 || i>1 && labels.(lblNames{j})(i)~=labels.(lblNames{j})(i-1)
                blkLabels{end+1}=mr.makeLabel('SET',lblNames{j},labels.(lblNames{j})(i));
            end
        end
        
        if ~isempty(blkLabels)
            iB=b_adc(i);
            b=seq.getBlock(iB);
            e=mr.block2events(b);
            seq.setBlock(iB,e{:},blkLabels{:});
        end
    end
    warning(warnBcp);

    % convert usable fields of aux to sequence definitions
    fieldsToExport={'kSpaceCenterLine','kSpaceCenterPartition','kSpaceCenterSample','kSpacePhaseEncodingLines',...
        'PhaseResolution','ReadoutOversamplingFactor','SliceGap','SlicePositions','SliceThickness',...
        'TargetGriddedSamples','TrapezoidGriddingParameters','AccelerationFactorPE','AccelerationFactor3D',...
        'FirstFourierLine','FirstRefLine','FirstFourier3D','FirstRef3D'};
    for i=1:length(fieldsToExport)
        if isfield(aux,fieldsToExport{i})
            prevDef=seq.getDefinition(fieldsToExport{i});
            if ~isempty(prevDef)
                warning('Overwriting existing sequence definition %s = %s', fieldsToExport{i}, num2str(prevDef));
            end
            seq.setDefinition(fieldsToExport{i},aux.(fieldsToExport{i})); 
        end
    end
end

%% test/plot label settings

if ~opt.skipApply
    lbls=seq.evalLabels('evolution','adc');
else
    lbls=labels;
end
lbl_names=fieldnames(lbls);
figure; hold on;
for n=1:length(lbl_names)
    plot(lbls.(lbl_names{n}));
end
legend(lbl_names(:));
title('evolution of labels/counters/flags');
xlabel('adc number');
