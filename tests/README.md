this is the directory containing test routines. you can either use them 
in Matlab or also run them from the command line as follows: 

1: change into the rood directory of the pulseq project (the one 
   containing 'doc', 'matlab', 'src' and 'tests' directories, amongst 
   others)
2: run the following command: 
   ```
   matlab -nodisplay -nosplash -r  "addpath(genpath('matlab')), results = runtests('tests'), exit(0);"
   ```
