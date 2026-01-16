function addCustomLabel(new_lbl)
% addCustomLabel(new_lbl)
% registers a new custom data label

if ~ischar(new_lbl)
    error('addCustomLabel: new label should be a character sctring');
end

supported_labels=mr.getSupportedLabels();
if any(ismember(supported_labels,new_lbl))
   warning('addCustomLabel: label %s is already known',new_lbl); 
end

supported_labels{end+1}=new_lbl;
mr.aux.globalVars('set','SupportedLabels',supported_labels);

end
