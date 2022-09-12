% this is an experimentaal high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping

% Set system limits
sys = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime',20e-6, 'B0', 2.89 ... % this is Siemens' 3T
);  

seq=mr.Sequence(sys);      % Create a new sequence object

%% prepare test objects
% pns is induced by the ramps, so we use long gradients to isolate the
% effects of the ramps
gpt=10e-3;
delay=30e-3;
rt_min=sys.gradRasterTime;
rt_test=floor(sys.maxGrad/sys.maxSlew/sys.gradRasterTime)*sys.gradRasterTime;
ga_min=sys.maxSlew*rt_min;
ga_test=sys.maxSlew*rt_test;

gx_min=mr.makeTrapezoid('x',sys,'amplitude',ga_min,'riseTime',rt_min,'fallTime',2*rt_min,'flatTime',gpt);
gy_min=gx_min;
gy_min.channel='y';
gz_min=gx_min;
gz_min.channel='z';

gx_test=mr.makeTrapezoid('x',sys,'amplitude',ga_test,'riseTime',rt_test,'fallTime',2*rt_test,'flatTime',gpt);
gy_test=gx_test;
gy_test.channel='y';
gz_test=gx_test;
gz_test.channel='z';

g_min={gx_min,gy_min,gz_min};
g_test={gx_test,gy_test,gz_test};

% dummy FID sequence
% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',0.1e-3, 'system', sys);
% Define delays and ADC events
adc = mr.makeAdc(512,'Duration',6.4e-3, 'system', sys);


%% Define sequence blocks
seq.addBlock(mr.makeDelay(delay));
for a=1:3
    seq.addBlock(g_min{a});
    seq.addBlock(mr.makeDelay(delay));
    seq.addBlock(g_test{a});
    seq.addBlock(mr.makeDelay(delay));
    for b=(a+1):3
        seq.addBlock(g_min{a},g_min{b});
        seq.addBlock(mr.makeDelay(delay));
        seq.addBlock(g_test{a},g_test{b});
        seq.addBlock(mr.makeDelay(delay));
    end
end
seq.addBlock(g_min);
seq.addBlock(mr.makeDelay(delay));
seq.addBlock(g_test);
seq.addBlock(mr.makeDelay(delay));
seq.addBlock(g_min);
seq.addBlock(rf);
seq.addBlock(adc);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% do some visualizations

seq.plot();             % Plot all sequence waveforms

%seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 25e-3]); %detailed view

%% 'install' to the IDEA simulator

seq.write('idea/external.seq');

%% PNS calc

%[pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
%[pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_GC98SQ.asc'); % aera-xq
[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc'); % TERRA-XR

%% load simulation results 

%[sll,~,~,vscale]=dsv_read('idea/dsv/prisma_pulseq_SLL.dsv');
%[sll,~,~,vscale]=dsv_read('idea/dsv/aera_pulseq_SLL.dsv');
[sll,~,~,vscale]=dsv_read('idea/dsv/terra_pulseq_SLL.dsv');
sll=cumsum(sll/vscale);

%% plot
figure;plot(sll(104:end)); % why 104? good question
hold on
plot(tpns*1e5-0.5,pns_n);
title('comparing internal and IDEA PNS predictions');

%% manual time alignment to calculate relative differences
ssl_s=104+tpns(1)*1e5-1.5;
ssl_e=ssl_s+length(pns_n)-1;
%figure;plot(sll(ssl_s:ssl_e));hold on; plot(pns_n);
figure;plot((sll(ssl_s:ssl_e)-pns_n)./pns_n*100);
title('relative difference in %');
