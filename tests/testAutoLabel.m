function tests = testAutoLabel
  if exist('functiontests')
    tests = functiontests(localfunctions);
  else
    lf=localfunctions();
    for i=1:length(lf)
      f=lf{i};
      n=func2str(f);
      if length(n)>3 && strcmp(n(1:4),'test')
        f();
        fprintf('Test function %s completed successfully\n', n);
      end
    end
  end
end

function test_detect_labels(testCase)
    % Build a simple GRE-like sequence with several ADC readouts
    seq = buildTestSeq(6);

    % Detect labels without applying them
    [labels, aux] = seq.autoLabel('skipApply', true);
    close all;

    if mr.aux.isOctave()
      assert(isstruct(labels));
    else
      testCase.verifyTrue(isstruct(labels));
    end

    % count ADCs in sequence
    nADCs = 0;
    for i = 1:length(seq.blockEvents)
        b = seq.getBlock(i);
        if isfield(b,'adc') && ~isempty(b.adc)
            nADCs = nADCs + 1;
        end
    end

    fld = fieldnames(labels);

    if mr.aux.isOctave()
      assert(~isempty(fld), 'autoLabel returned no label fields');
    else
      testCase.verifyNotEmpty(fld, 'autoLabel returned no label fields');
    end

    % Each label vector must match number of ADCs and contain non-negative integers
    for i = 1:length(fld)
        v = labels.(fld{i});
        if mr.aux.isOctave()
          assert(length(v)==nADCs, sprintf('Label %s has wrong length', fld{i}));
          assert(min(v)>=0, sprintf('Label %s contains negative values', fld{i}));
          assert(all(abs(v - round(v)) < 1e-10), sprintf('Label %s contains non-integer values', fld{i}));
        else
          testCase.verifyEqual(length(v), nADCs, sprintf('Label %s has wrong length', fld{i}));
          testCase.verifyGreaterThanOrEqual(min(v), 0, sprintf('Label %s contains negative values', fld{i}));
          testCase.verifyTrue(all(abs(v - round(v)) < 1e-10), sprintf('Label %s contains non-integer values', fld{i}));
        end
    end

    % aux should contain at least center sample info
    if mr.aux.isOctave()
      assert(isfield(aux, 'kSpaceCenterSample') || isempty(fieldnames(aux))==0);
    else
      testCase.verifyTrue(isfield(aux, 'kSpaceCenterSample') || isempty(fieldnames(aux))==0);
    end
end

function test_apply_labels(testCase)
    % Build sequence and detect labels
    seq = buildTestSeq(8);
    [labels, ~] = seq.autoLabel('skipApply', true);
    close all;

    % Now apply labels to the sequence
    seq.autoLabel();
    close all;

    % Evaluate labels from sequence after apply
    lbls = seq.evalLabels('evolution', 'adc');

    f = fieldnames(labels);
    for i = 1:length(f)
      if mr.aux.isOctave()
        assert(isfield(lbls, f{i}), sprintf('Applied labels missing field %s', f{i}));
        assert(all(lbls.(f{i})==labels.(f{i})), sprintf('Applied label %s does not match detected values', f{i}));
      else
        testCase.verifyTrue(isfield(lbls, f{i}), sprintf('Applied labels missing field %s', f{i}));
        testCase.verifyEqual(lbls.(f{i}), labels.(f{i}), sprintf('Applied label %s does not match detected values', f{i}));
      end
    end
end

function test_reflect_reorder_3d(testCase)
    % Build a simple 3D Cartesian-like sequence
    seq = buildTestSeq3D(4, 3);

    % Baseline labels
    labels0 = seq.autoLabel('skipApply', true);
    close all;

    if mr.aux.isOctave()
      assert(isfield(labels0, 'LIN'));
      assert(isfield(labels0, 'PAR'));
    else
      testCase.verifyTrue(isfield(labels0, 'LIN'));
      testCase.verifyTrue(isfield(labels0, 'PAR'));
    end

    % Reorder should swap encoding axes for this sequence
    labelsReorder = seq.autoLabel('skipApply', true, 'reorder', [1 3 2]);
    close all;
    if mr.aux.isOctave()
      assert(all(labelsReorder.LIN==labels0.PAR), ...
        'reorder did not swap LIN/PAR as expected');
      assert(all(labelsReorder.PAR==labels0.LIN), ...
        'reorder did not swap PAR/LIN as expected');
    else
        testCase.verifyEqual(labelsReorder.LIN, labels0.PAR, ...
        'reorder did not swap LIN/PAR as expected');
      testCase.verifyEqual(labelsReorder.PAR, labels0.LIN, ...
        'reorder did not swap PAR/LIN as expected');
    end

    % Reflect should reverse counters along reflected axes
    labelsReflect = seq.autoLabel('skipApply', true, 'reflect', [2 3]);
    close all;
    if mr.aux.isOctave()
      assert(all(labelsReflect.LIN==max(labels0.LIN) - labels0.LIN), ...
        'reflect did not reverse LIN as expected');
      assert(all(labelsReflect.PAR==max(labels0.PAR) - labels0.PAR), ...
        'reflect did not reverse PAR as expected');
    else
      testCase.verifyEqual(labelsReflect.LIN, max(labels0.LIN) - labels0.LIN, ...
        'reflect did not reverse LIN as expected');
      testCase.verifyEqual(labelsReflect.PAR, max(labels0.PAR) - labels0.PAR, ...
        'reflect did not reverse PAR as expected');
    end
end

function seq = buildTestSeq(nLines)
    % Helper to build a simple GRE-like sequence with nLines ADC readouts
    seq = mr.Sequence();
    % small slice-selective excitation (provide slice gradient)
    [rf0, gz0] = mr.makeSincPulse(pi/8, 'duration', 2e-3, 'SliceThickness', 5e-3, 'use', 'excitation');
    seq.addBlock(rf0, gz0);
    for i = 1:nLines
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000));
        seq.addBlock(mr.makeTrapezoid('y', 'area', -500 + (i-1) * 100));
        seq.addBlock(mr.makeTrapezoid('x', 'area', -500));
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000, 'duration', 10e-3), mr.makeAdc(64, 'duration', 10e-3));
    end
end

function seq = buildTestSeq3D(nLines, nPartitions)
    % Helper to build a simple 3D Cartesian-like sequence
    seq = mr.Sequence();

    [rf0, gz0] = mr.makeSincPulse(pi/8, 'duration', 1e-3, 'SliceThickness', 100e-3, 'use', 'excitation');
    
    for ip = 1:nPartitions
        for il = 1:nLines
            kyArea = -300 + (il-1) * 200;
            kzArea = -200 + (ip-1) * 200;

            seq.addBlock(rf0, gz0);
            seq.addBlock(mr.makeTrapezoid('x', 'area', 400));
            seq.addBlock(mr.makeTrapezoid('y', 'area', kyArea), mr.makeTrapezoid('z', 'area', kzArea));
            seq.addBlock(mr.makeTrapezoid('x', 'area', 800, 'duration', 8e-3), mr.makeAdc(64, 'duration', 8e-3));
            seq.addBlock(mr.makeTrapezoid('x', 'area', -1200));
        end
    end
end
