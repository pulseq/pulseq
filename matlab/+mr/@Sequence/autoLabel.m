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
%   autoLabel(...,'useLabels',labels_struct) Skip label evolution
%   detection and only apply externally-calculated labels to the sequence.
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
    parser.addParamValue('skipApply',false,@(x)(islogical(x) && isscalar(x)));
end
parse(parser,varargin{:});
opt = parser.Results;

% if isempty(opt.useLabels)
%     labels=struct();
% else
%     labels=opt.useLabels;
% end

if ~isfinite(opt.blockRange(2))
    opt.blockRange(2)=length(seq.blockEvents);
end

aux = struct();

if isempty(opt.useLabels)

    %% preparing calculations
    [ktraj_adc, t_adc, ~, ~, t_excitation, ~, slicepos, t_slicepos, gw_pp] = seq.calculateKspacePP('blockRange',opt.blockRange);
    blockStartTimes=cumsum([0 seq.blockDurations]);
    blockStartTimes=blockStartTimes(opt.blockRange(1):opt.blockRange(2));
    
    %% slice positions
    sliceGrads=zeros(size(slicepos));
    for i=1:3
        sliceGrads(i,:)=ppval(gw_pp{i},t_slicepos);
    end
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
    % useful results:
    %   uniqueSlicePositions : speaks for itself, sorted in the order of occurence
    %   sliceCountersAcquisitionOrder : slice cointers in the acquisition order
    %   b_usp: block indices containg RF pulses corresponding to the first
    %          occurence of each of the unique slice positions
    %   sortedSlicePositions : slice positions in the accending order
    %   sliceCountersSorted : slice counters corresponding to the sorted slice order
    
    %% index ADCs
    adcLengths=[];
    for iB=opt.blockRange(1):opt.blockRange(2)
        block = seq.getBlock(iB);
        if ~isempty(block.adc)
            adcLengths(end+1)=block.adc.numSamples;
        end
    end
    adcStartCounts=cumsum([1, adcLengths(1:end-1)]);
    t_adcStarts=t_adc(adcStartCounts);
    b_adc=zeros(size(t_adcStarts));
    for i=1:numel(t_adcStarts)
        b_adc(i)=find(blockStartTimes<t_adcStarts(i),1,'last');
    end
    firstNonNoiseAdc=find(t_adcStarts>t_excitation(1),1,'first');
    firstNonNoiseAdcSampe=find(t_adc>t_excitation(1),1,'first');
    sliceCountersAdc=zeros(size(t_adcStarts)-firstNonNoiseAdc+1);
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
    nReadouts = numel(t_adcStarts)-firstNonNoiseAdc+1;
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
    sampler=zeros(kindex_end');
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
        repeat(i)=sampler(ind);
        sampler(ind)=repeat(i)+1;
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
else
    labels=opt.useLabels;
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
