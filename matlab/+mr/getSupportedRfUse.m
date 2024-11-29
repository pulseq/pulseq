function [supported_rf_use, short_rf_use] = getSupportedRfUse()
% auxilary function

supported_rf_use={'excitation','refocusing','inversion','saturation','preparation','other','undefined'};
if nargout>1
    short_rf_use=cell2mat(cellfun (@(x) x(1),supported_rf_use,'un',0));
end

end
