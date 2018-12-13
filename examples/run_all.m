% run all example sequences
%
% Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

% avoid cwd madness
od = cd(fileparts(which(mfilename)));
addpath(genpath('../matlab'));

% list of sequences to execute
fid
gre
epi_rs

% just in case the change of the cwd persists in the parent shell
cd(od);