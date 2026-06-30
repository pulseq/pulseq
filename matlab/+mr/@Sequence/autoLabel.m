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
%   specified) prior to generating label ranges. Affects both Fourier
%   encoding dimensions (redout, phase and partition encoding) and the
%   slice ordering. If used in combination with 'reorder', it is performed
%   first. 
%
%   autoLabel(...,'reorder',[2,1]) Reorder axes of k-space trajectories
%   (in this example by swapping directions 1 and 2 prior to generating
%   label ranges (reordering vector can contain 2 or 3 entries). If used in
%   combination with 'reflect', it is performed last.
% 
%   autoLabel(...,'mirrorFourier',true) Mirrors all Fourier encoding
%   directions, e.g. if ifft() is used instead of fft() by the image
%   reconstruction algorithm. It essentially reverts all Fourier encoding
%   dimensions (readout, phase, partition encoding) but does not affect
%   slice ordering. Can be arbitrarily combined with 'reflect'.
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
%   From the experience, Siemens scanners need 'mirrorFourier' if the "XYZ
%   in TRA" coordinate mapping (default) mode is applied, which transforms
%   the sequence by the following rotation matrix [0 -1 0; -1 0 0; 0 0 -1]
%   prior to passing it to IDEA; On Siemens 'sortSlices'='descending' is
%   optional, otherwise the interpreter will change the slice indexes.
%
%   see evalLabels()

validSliceSorting={'acquisition','ascending','descending'};
persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'autoLabel';
    parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
    parser.addParamValue('useLabels',struct([]),@(x)(isempty(x) || isstruct(x)));
    parser.addParamValue('useAux',struct([]),@(x)(isempty(x) || isstruct(x)));    
    parser.addParamValue('skipApply',false,@(x)(islogical(x) && isscalar(x)));
    parser.addParamValue('mirrorFourier', false, @(x)(islogical(x) && isscalar(x)));
    parser.addParamValue('reflect', [], @(x)(isnumeric(x) && numel(x)<4));
    parser.addParamValue('reorder', [], @(x)(isnumeric(x) && numel(x)<4));
    parser.addParamValue('sortSlices','acquisition',@(x) any(validatestring(x,validSliceSorting)));
    parser.addParamValue('noPlots', false, @(x)(islogical(x) && isscalar(x)));
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
    raw_block = seq.getRawBlockContentIDs(iB); % this is much faster than seq.getBlock(iB)
    if ~isempty(raw_block.adc)
        libData = seq.adcLibrary.data(raw_block.adc).array;
        % from getBlock(): adc = cell2struct(num2cell(libData(1:end-1)), {'numSamples', 'dwell', 'delay', ...
        adcLengths(end+1)=libData(1); % libData(1) is adc.numSamples
        b_adc(end+1)=iB-opt.blockRange(1)+1;
    end
    %block = seq.getBlock(iB);
    %if ~isempty(block.adc)
    %    adcLengths(end+1)=block.adc.numSamples;
    %    b_adc(end+1)=iB-opt.blockRange(1)+1;
    %end
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
    
    if opt.mirrorFourier
        ktraj_adc=-ktraj_adc;
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
    [sortedSlicePositions, ~, sliceCountersSortedPositions] = unique(sliceOffsets);

    % aux.SliceThickness
    % aux.SliceGap
    % todo: more checks, e.g. if all slice thickeness and gaps are the same...
    GAs1=vecnorm(sliceGrads(:,1));
    if GAs1>0
        bSlice1=find(blockStartTimes<t_slicepos(1),1,'last');
        b=seq.getBlock(bSlice1);
        bw=mr.calcRfBandwidth(b.rf);
        aux.SliceThickness=bw/GAs1;
        if length(sortedSlicePositions)>1
            aux.SliceGap=diff(sortedSlicePositions(1:2))-aux.SliceThickness;
        end
    end
    if strcmp(opt.sortSlices,'descending')
        sortedSlicePositions=sortedSlicePositions(end:-1:1);
        sliceCountersSortedPositions = max(sliceCountersSortedPositions) + 1 - sliceCountersSortedPositions;
    end

    % useful results:
    %   uniqueSlicePositions : speaks for itself, sorted in the order of occurence
    %   sliceCountersAcquisitionOrder : slice cointers in the acquisition order
    %   b_usp: block indices containg RF pulses corresponding to the first
    %          occurence of each of the unique slice positions
    %   sortedSlicePositions : slice positions in the ascending or descending deorder
    %   sliceCountersSorted : slice counters corresponding to the sorted slice order
    
    %% index ADCs, part 2
    if ~strcmp(opt.sortSlices,'acquisition')
        uniqueSlicePositions=sortedSlicePositions;
    end
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
    tEcho = zeros(1,nReadouts);
    kEcho = zeros(3,nReadouts);
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
        kEcho(:,i)=ktraj_adc(:,c1+cEchoPos(i)-1);
        t_adcThisEcho=t_adc(c1+cEchoPos(i)-1);
        if vecnorm(kspaceCenterPoint)>eps % this is more or less copied from seq.testRepoert()
            % the actual echo might be between k-space samples, try to interpolate it
            i2check=[];
            % check if adc kspace trajectory has elements left and right to index_echo
            if cEchoPos(i) > 1
                i2check=c1+cEchoPos(i)-2;
            end
            if c1+cEchoPos(i)-1 < c2
                i2check(end+1)=c1+cEchoPos(i);
            end
            for a=1:numel(i2check)
                v_i_to_0=-kEcho(:,i);
                v_i_to_t=ktraj_adc(:,i2check(a))-kEcho(:,i);
                % project v_i_to_0 to v_o_to_t
                p_vit=v_i_to_0'*v_i_to_t/(vecnorm(v_i_to_t)^2);
                if p_vit>0
                    % we have forund a bracket for the echo and the proportionality coefficient is p_vit
                    t_adcThisEcho=t_adcThisEcho*(1-p_vit) + t_adc(i2check(a))*p_vit;
                    break;
                end
            end
        end

        i_sliseposThisEcho=find(t_slicepos<t_adcThisEcho,1,'last');
        tEcho(i)=t_adcThisEcho-t_slicepos(i_sliseposThisEcho);
        for j=1:3
            gradReadout(j,i) = ppval(gw_pp{j},t_adcThisEcho); % this call is relatively expensive for long sequences, we should go away from it       
        end
        if opt.mirrorFourier
            gradReadout(:,i)=-gradReadout(:,i);
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
        if i>1 && vecnorm(gradReadout(:,i)*signReadout(i)-gradReadout(:,1)*signReadout(1))>1e-4 % why such threshold? 
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

    % detect non-constant sampling
    if any(abs(diff(kCentralReadoutProjection,2)/dkR*numel(kCentralReadoutProjection))>0.1)
        % see if we can extraxt trapezoid resampling parameters
        bCentralRO=find(blockStartTimes<t_adcStarts(firstNonNoiseAdcSampe+cCentralReadout-1),1,'last');
        b=seq.getBlock(bCentralRO);
        if ~isempty(opt.reorder) && opt.reorder(1)~=1
            warning('EPI readout gridding parameter detection code is incompatible with the reorder option');
        else
            if isfield(b,'gx') && ~isempty(b.gx) && strcmp(b.gx.type,'trap') && isfield(b,'adc') && ~isempty(b.adc)
                aux.TrapezoidGriddingParameters=[b.gx.riseTime, b.gx.flatTime, b.gx.fallTime, b.adc.delay-b.gx.delay b.adc.numSamples*b.adc.dwell];
                % TODO: figure out a better way to find TargetGriddedSamples... In Siemens sequences it is often adc.numsamples
                aux.TargetGriddedSamples=b.adc.numSamples;
            end
        end
        % TODO: figure out a better way to find TargetGriddedSamples... In Siemens sequences it is often adc.numsamples
    end

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
        if ~opt.noPlots
            figure; plot(aux.kReadout(1,:),'.-'); hold on; plot(aux.kReadout(2,:),'.-');plot(aux.kReadout(2,end:-1:1),'x'); title('bipolar readout k-space alignment'); 
            legend('forward','reverse','rev.refl'); % if crosses overlap the line then REV is symmetric
            xlabel('sample number (acq. order)');
            ylabel('k-space, 1/m');
        end
    else
        aux.kReadout = kCentralReadoutProjection;
    end

    % useful results: 
    %   isCartesianReadout : is this trajectory sufficiently Cartesian-like
    %   gradReadout : readout gradients as 3D vectors
    %   signReadout : sign of the dominant readout component
    %   kCentralReadout : the most central readout
    %   cCentralReadoutCenter : center of the above
    %   cEchoPos : position of the echo within each ADC vector for each readout
    %   tEcho : echo time for each readout
    %   kEcho : k-space location of each echo

    %% detect navigators
    isNavigator=false(1,nReadouts);
    isNavigatorCandidate = vecnorm(kEcho-kspaceCenterPoint) < 1e-4; % why this threshold?
    orderedReadoutIndicator=diff(kEcho-kspaceCenterPoint,1,2);
    i=find(max(abs(orderedReadoutIndicator'))>1e-4);
    orderedReadoutIndicator=orderedReadoutIndicator(i,:); % discard empty dimensions
    %TODO: finish this detection some day... e.g. if all(isNavigatorCandidate(1:3)) && abs(orderedReadoutIndicator(1)-orderedReadoutIndicator(2))<1e-4 && max(abs(diff(orderedReadoutIndicator(4:x))))<1e-4
    if nReadouts>=16 && all(isNavigatorCandidate(1:3)) && abs(orderedReadoutIndicator(1)-orderedReadoutIndicator(2))<1e-4 && ...
            max(abs(diff(orderedReadoutIndicator(4:16))))<1e-4
        aux.epiWithThreeEchoNavigator=true;
        % this is a real hack to close those singular k=0 lines in the middle of readouts
        isNavigator = (isNavigatorCandidate + circshift(isNavigatorCandidate,1) + circshift(isNavigatorCandidate,-1))>1.5;
        % TODO: set SEG (010 10101...) and AVG (001 00000...)
    end

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
        if (~isNavigator(i))
            % ignore navigators
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
    end
    % if (max(repeat(:))>0)
    %     kindex=[kindex;(repeat+1)];
    %     kindex_mat=[kindex_mat;(repeat+1)];
    %     kindex_end=max(kindex_mat,[],2);
    % end
    %figure; plot(kindex(1,:),kindex(2,:),'.-');
    %figure; plot(kindex_mat(1,:),'.-');

    % see if some repetitions are actully echoes/contrasts
    nRep=max(repeat)+1;
    skipECO=false;
    if nRep>1
        TE=zeros(1, nRep);
        for i=1:nRep
            try
                if isvector(kindex)
                    TE(i)=tEcho(kindex==0 & sliceCountersAdc==1 & repeat==(i-1)); % we assume all slices have the same TEs and only check slice 1
                else
                    TE(i)=tEcho(all(kindex==0) & sliceCountersAdc==1 & repeat==(i-1)); % we assume all slices have the same TEs and only check slice 1
                end
            catch
                warning('unclear sequence structure, skipping TE & ECO counter detection');
                skipECO=true;
                break;
            end
        end
        if ~skipECO
            [TE_sorted, TE_order]=sort(TE);
            TE_cluster=cumsum([1, diff(TE_sorted)>10e-6]);
            unique_TE=zeros(1,max(TE_cluster));
            for i=1:length(unique_TE)
                unique_TE(i)=mean(TE_sorted(TE_cluster==i));
                TE(TE_order(TE_cluster==i))=unique_TE(i);
            end
            aux.TE=unique_TE; % maybe we should fill this also for single-TE sequences?
            echo=zeros(1,size(ktraj_adc,2));
            echo_rep=zeros(1,nRep);
            for i=1:nRep
                cecho=find(abs(unique_TE-TE(i))<=1e-6);
                echo(repeat==(i-1))=cecho;
                repeat(repeat==(i-1))=sum(echo_rep==cecho);
                echo_rep(i)=cecho;
            end
        end
    end

    %% create labels struct
    labels=struct();
    % NOISE
    % SLC
    % REV
    % LIN,
    % PAR
    % ECO
    % REP 
    
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
    if exist('echo','var') && max(echo(:))>1 
        if firstNonNoiseAdc>1
            labels.ECO=zeros(1,nADCs);
            labels.ECO(firstNonNoiseAdc:end)=echo(:).'-1;
        else
            labels.ECO=echo(:).'-1;
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
    if any(isNavigator)
        labels.NAV=isNavigator;
    end

    if isfield(labels,'LIN'), aux.kSpaceCenterLine = labels.LIN(cCentralReadout); end
    if isfield(labels,'PAR'), aux.kSpaceCenterPartition = labels.PAR(cCentralReadout); end
    aux.kSpaceCenterSample = cCentralReadoutCenter-1;
    if length(uniqueSlicePositions)>1
        aux.SlicePositions = uniqueSlicePositions;
    end
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
if opt.noPlots
    return;
end

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
