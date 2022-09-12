function [ report ] = testReport( obj )
%testReport Analyze the sequence and return a text report
%   Currently no parameters are required. In future versions it may be
%   possible to (de)select some tests
%
% maxim.zaitsev@uniklinik-freiburg.de

% find the RF pulses and list flip angles
flipAnglesDeg=[];
for k=obj.rfLibrary.keys
    libData=obj.rfLibrary.data(k).array;
    if length(obj.rfLibrary.type)>=k
        rf=obj.rfFromLibData(libData,obj.rfLibrary.type(k));
    else
        rf=obj.rfFromLibData(libData);
    end
    flipAnglesDeg=[flipAnglesDeg abs(sum(rf.signal(1:end-1).*(rf.t(2:end)-rf.t(1:end-1))))*360]; %we use rfex.t(1) in place of opt.system.rfRasterTime
end
flipAnglesDeg=unique(flipAnglesDeg);

% calculate TE and TR

[duration, numBlocks, eventCount]=obj.duration();

%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = obj.calculateKspace();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = obj.calculateKspacePP();
[ktraj_adc, t_adc, ~, ~, t_excitation, ~] = obj.calculateKspacePP();

% trajectory calculation will fail for spin-echoes if seq is loaded from a 
% file for the current file format revision (1.2.0) because we do not store 
% the use of the RF pulses. Read function has an option 'detectRFuse' which
% may help...

%        
kabs_adc=sum(ktraj_adc.^2,1).^0.5;
[kabs_echo, index_echo]=min(kabs_adc);
t_echo=t_adc(index_echo); % just a first estimate, see if we can improve it
if kabs_echo>eps
    i2check=[];
    % check if adc kspace trajectory has elements left and right to index_echo
    if index_echo > 1
        i2check=[i2check (index_echo-1)];
    end
    if index_echo < length(kabs_adc)
        i2check=[i2check (index_echo+1)];
    end
    for a=1:numel(i2check)
        v_i_to_0=-ktraj_adc(:,index_echo);
        v_i_to_t=ktraj_adc(:,i2check(a))-ktraj_adc(:,index_echo);
        % project v_i_to_0 to v_o_to_t
        p_vit=v_i_to_0'*v_i_to_t/(vecnorm(v_i_to_t)^2);
        if p_vit>0
            % we have forund a bracket for the echo and the proportionality
            % coefficient is p_vit
            t_echo=t_adc(index_echo)*(1-p_vit) + t_adc(i2check(a))*p_vit;
            break;
        end
    end
end

if ~isempty(t_excitation)
    t_ex_tmp=t_excitation(t_excitation<t_echo);
    TE=t_echo-t_ex_tmp(end);
    % TODO detect multiple TEs
else
    TE=NaN;
end

if (length(t_excitation)<2)
    TR=duration; % best estimate for now
else
    t_ex_tmp1=t_excitation(t_excitation>t_echo);
    if isempty(t_ex_tmp1)
        TR=t_ex_tmp(end)-t_ex_tmp(end-1);
    else
        TR=t_ex_tmp1(1)-t_ex_tmp(end);
    end
    % TODO check frequency offset to detect multiple slices
end

% check sequence dimensionality and spatial resolution
k_extent=max(abs(ktraj_adc),[],2);
k_scale=max(k_extent);
if (k_scale~=0)
    k_bins=4e6; % this defines our ability to separate k-space samples. 
                % lower values give us imunity to rounding errors in k-space calculations
                % current code below (2nd pass) however merges neighboring cells (+-1)
    k_threshold=k_scale/k_bins;

    % detect unused dimensions and delete them
    if any(k_extent<k_threshold)
        ktraj_adc(k_extent<k_threshold,:)=[]; % delete rows
        k_extent(k_extent<k_threshold)=[];
    end

    % bin the k-space trajectory to detect repetitions / slices
    k_len=size(ktraj_adc,2);
    k_repeat=zeros(1,k_len);
    k_storage=zeros(1,k_len);
    k_storage_next=1;
    % containers.Map only supports string as a key... (in older matlabs)
    %kmap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    kmap = java.util.HashMap;
    for i=1:k_len
        key_string = sprintf('%d ', int32(k_bins+round(ktraj_adc(:,i)/k_threshold))); 
        % containers.Map does not have a proper find function so we use direct
        % access and catch the possible error
        k_storage_ind = kmap.get(key_string);
        if isempty(k_storage_ind)
            k_storage_ind=k_storage_next;
            kmap.put(key_string,k_storage_ind);
            k_storage_next=k_storage_next+1;
        end
        k_storage(k_storage_ind)=k_storage(k_storage_ind)+1;
        k_repeat(i) = k_storage(k_storage_ind);
    end
    % at this point k_storage(1:(k_storage_next-1)) is our visit frequency map
    Repeats_max=max(k_storage(1:(k_storage_next-1)));
    Repeats_min=min(k_storage(1:(k_storage_next-1)));
    Repeats_median=median(k_storage(1:(k_storage_next-1))); 
    Repeats_unique=unique(k_storage(1:(k_storage_next-1)));
    Counts_unique=zeros(size(Repeats_unique));
    for i=1:numel(Repeats_unique)
        Counts_unique(i)=sum(Repeats_unique(i)==k_storage(1:(k_storage_next-1)));
    end

    ktraj_rep1=ktraj_adc(:,k_repeat==1);
    % TODO: think of something clever, e.g. detecting maximum delta-k
    % if length(k_extent)==2
    %     dt = delaunayTriangulation(ktraj_rep1');
    %     k = convexHull(dt);
    %     figure; plot(dt.Points(:,1),dt.Points(:,2), '.', 'markersize',10); hold on;
    %     plot(dt.Points(k,1),dt.Points(k,2), 'r'); hold off;
    %     %[V,R] = voronoiDiagram(dt);
    %     figure; voronoi(dt);
    % end
    % try to detect k-space lines or columns.
    k_counters=zeros(size(ktraj_rep1));
    dims=size(ktraj_rep1,1);
    %ordering=cell(1,dims);
    kmap = java.util.HashMap;
    for j=1:dims
        %kmap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
        k_storage=zeros(1,k_len);
        k_storage_next=1; 
        kmap.clear();
        for i=1:size(ktraj_rep1,2)
            key=int32(round(ktraj_rep1(j,i)/k_threshold));
            k_storage_ind = kmap.get(key);
            if isempty(k_storage_ind) % attempt to account for rounding errors
                k_storage_ind = kmap.get(key+1);
            end
            if isempty(k_storage_ind) % attempt to account for rounding errors
                k_storage_ind = kmap.get(key-1);
            end
            if isempty(k_storage_ind)
                k_storage_ind=k_storage_next;
                kmap.put(key,k_storage_ind);
                k_storage_next=k_storage_next+1;
                k_storage(k_storage_ind)=ktraj_rep1(j,i);
                %fprintf('%d:%d(%g) ',k_storage_ind,key,ktraj_rep1(j,i));
            end
            %assert(k_storage_ind==k_storage(k_storage_ind));
            k_counters(j,i) = k_storage_ind;
        end
        %ordering{j}=cell2mat(kmap.values);
        %fprintf('\n');
    end
    unique_kpositions=max(k_counters,[],2);
    isCartesian=(prod(unique_kpositions)==size(ktraj_rep1,2));
else
    unique_kpositions=1;
end
% check gradient amplitudes and slew rates

% gradient waveform
gw_data=obj.waveforms_and_times(); % FIXME: avoid this second call for generating gradient shapes (1st one was inside of the k-space calculation routine)
gws=cell(size(gw_data));
ga=zeros(length(gw_data),1);
gs=zeros(length(gw_data),1);
% to calculate max absolute gradients and slew rates we have to play
% tricks... we namely have to interpolate the data to the common time axis
dim1ind = @(x, n) x(n,:);
common_time=unique(dim1ind([gw_data{:}],1));
gw_ct=zeros(length(gw_data),length(common_time));
gs_ct=zeros(length(gw_data),length(common_time)-1);
for gc=1:length(gw_data)
    if size(gw_data{gc},2)>0 
        gws{gc}=(gw_data{gc}(2,2:end)-gw_data{gc}(2,1:end-1))./(gw_data{gc}(1,2:end)-gw_data{gc}(1,1:end-1)); % slew
        % interpolate to the common time
        gw_ct(gc,:)=interp1(gw_data{gc}(1,:),gw_data{gc}(2,:),common_time,'linear',0);
        gs_ct(gc,:)=(gw_ct(gc,2:end)-gw_ct(gc,1:end-1))./(common_time(2:end)-common_time(1:end-1));
        % max grad/slew per channel
        ga(gc)=max(abs(gw_data{gc}(2,:)));
        gs(gc)=max(abs(gws{gc}));
        % TODO: calculate grad RMS values (this is an interesting task in the piece-wise-linear domain)
    end
end

%figure; plot(common_time, gw_ct');
%figure; plot(common_time, sum(gw_ct.^2,1).^0.5);
%figure; plot(0.5*(common_time(1:end-1)+common_time(2:end)), sum(gs_ct.^2,1).^0.5);

% max absolute value grad/slew -- check for a worst case upon rotation
ga_abs=max(sum(gw_ct.^2,1).^0.5);
gs_abs=max(sum(gs_ct.^2,1).^0.5);


% check timing of blocks and delays (raster alignment)
[timing_ok, timing_error_report] = obj.checkTiming();

report = { sprintf('Number of blocks: %d\n',numBlocks),...
           [ sprintf('Number of events:\n'),...
             sprintf('   RF:   %6d\n',eventCount(2)),...
             sprintf('   Gx:   %6d\n',eventCount(3)),...
             sprintf('   Gy:   %6d\n',eventCount(4)),...
             sprintf('   Gz:   %6d\n',eventCount(5)),...
             sprintf('   ADC:  %6d\n',eventCount(6)),...
             sprintf('   Delay:%6d\n',eventCount(1)) ],...
           [ sprintf('Sequence duration: %.6fs\n',duration),...
             sprintf('TE: %.6fs\n',TE),...
             sprintf('TR: %.6fs\n',TR) ],...
             sprintf('Flip angle: %.02f°\n', flipAnglesDeg),...
             sprintf('Unique k-space positions (a.k.a. columns, rows, etc): %d\n', unique_kpositions)};
if any(unique_kpositions>1)
    report = { report{:},...
           [ sprintf('Dimensions: %d\n', length(k_extent)),...
             sprintf('   Spatial resolution: %.02f mm\n', 0.5./k_extent*1e3) ],...
             sprintf('Repetitions/slices/contrasts: %.d  range: [%.d %.d]\n', Repeats_median, Repeats_min, Repeats_max) };
    report = { report{:},...
        sprintf('   %d k-space position(s) repeated %d times\n', [Counts_unique;Repeats_unique])};

    if isCartesian
       report = { report{:}, sprintf('Cartesian encoding trajectory detected\n') };
    else
       report = { report{:}, sprintf('Non-Cartesian/irregular encoding trajectory detected (e.g. EPI, spiral, radial, etc)\n') };
    end
end
if (timing_ok)
    report = { report{:}, sprintf('Block timing check passed successfully\n') };
else
    report = { report{:}, [ sprintf('Block timing check failed! Error listing follows:\n'),...
                            sprintf([timing_error_report{:}]) ] };
end
report = { report{:},...
    sprintf('Max. Gradient: %.0f Hz/m == %.02f mT/m\n', [ga mr.convert(ga,'Hz/m','mT/m')]'),...
    sprintf('Max. Slew Rate: %g Hz/m/s == %.02f T/m/s\n', [gs mr.convert(gs,'Hz/m/s','T/m/s')]') };
report = { report{:},...
    sprintf('Max. Absolute Gradient: %.0f Hz/m == %.02f mT/m\n', [ga_abs mr.convert(ga_abs,'Hz/m','mT/m')]'),...
    sprintf('Max. Absolute Slew Rate: %g Hz/m/s == %.02f T/m/s\n', [gs_abs mr.convert(gs_abs,'Hz/m/s','T/m/s')]') };

end

