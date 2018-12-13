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
    rf=obj.rfFromLibData(libData);
    flipAnglesDeg=[flipAnglesDeg abs(sum(rf.signal))*rf.t(1)*360]; %we use rfex.t(1) in place of opt.system.rfRasterTime
end
flipAnglesDeg=unique(flipAnglesDeg);

% calculate TE and TR

[duration, numBlocks, eventCount]=obj.duration();

[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = obj.calculateKspace();

% trajectory calculation will fail for spin-echoes if seq is loaded from a 
% file for the current file format revision (1.2.0) because we do not store 
% the use of the RF pulses. Read function has an option 'detectRFuse' which
% may help...

%        
kabs_adc=sum(ktraj_adc.^2,1).^0.5;
[kabs_echo, index_echo]=min(kabs_adc);
t_echo=t_adc(index_echo);
t_ex_tmp=t_excitation(t_excitation<t_echo);
TE=t_echo-t_ex_tmp(end);
% TODO detect multiple TEs

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
    k_bins=4e6;
    k_threshold=k_scale/k_bins;

    % detect unused dimensions and delete them
    if any(k_extent<k_threshold)
        ktraj_adc(k_extent<k_threshold,:)=[]; % delete rows
        k_extent(k_extent<k_threshold)=[];
    end

    % bin the k-space trajectory to detect repetitions / slices
    k_len=size(ktraj_adc,2);
    k_repeat=zeros(1,k_len);
    % containers.Map only supports string as a key...
    kmap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    for i=1:k_len
        key_string = sprintf('%d ', k_bins+round(ktraj_adc(:,i)/k_threshold)); 
        % containers.Map does not have a proper find function so we use direct
        % access and catch the possible error
        try
            k_repeat(i) = kmap(key_string)+1;
        catch 
            k_repeat(i) = 1;
        end
        kmap(key_string) = k_repeat(i);
    end
    Repeats=max(k_repeat);

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
    ordering=cell(1,dims);
    for j=1:dims
        c=1;
        kmap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
        for i=1:size(ktraj_rep1,2)
            %key_string = sprintf('%d %d', dims(dims~=j), round(ktraj_rep1(dims~=j,i)/k_threshold));
            key=round(ktraj_rep1(j,i)/k_threshold);
            try
                k_counters(j,i) = kmap(key);
            catch
                k_counters(j,i) = c;
                kmap(key) = c;
                c=c+1;
            end
        end
        ordering{j}=cell2mat(kmap.values);
    end
    unique_kpositions=max(k_counters,[],2);
    isCartesian=(prod(unique_kpositions)==size(ktraj_rep1,2));
else
    unique_kpositions=1;
end
% check gradient amplitudes and slew rates

gw=obj.gradient_waveforms();
ga=max(abs(gw),[],2);
gs=max(abs(gw(:,2:end)-gw(:,1:end-1)),[],2)./obj.sys.gradRasterTime;

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
if (unique_kpositions>1)
    report = { report{:},...
           [ sprintf('Dimensions: %d\n', length(k_extent)),...
             sprintf('   Spatial resolution: %.02f mm\n', 0.5./k_extent*1e3) ],...
             sprintf('Repetitions/slices/contrasts: %.d\n', Repeats) };
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

end

