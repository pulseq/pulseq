% example script to plot gradient frequency spectrum
% expects: 
%  "seq" object to be already populated 
% furthermore, on Siemens you need the *.asc file for your gradient system

ascName=[]; % this disables the display of the system's resonance frequences
%ascName='idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'; % 3T prisma
%ascName='idea/asc/MP_GradSys_P034_X60.asc'; % 3T cima.X
%ascName='idea/asc/MP_GradSys_K2368_2250V_1250A_XT_GC25.asc'; % 3T vida

seq.gradSpectrum(ascName);
