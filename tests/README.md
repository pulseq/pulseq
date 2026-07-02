This is the directory containing test routines. You can either use them 
in Matlab or also run them from the command line as follows: 

1: change into the root directory of the pulseq project (the one 
   containing 'doc', 'matlab', 'src' and 'tests' directories, amongst 
   others)
2: run the following command: 
   ```
   matlab -nodisplay -nosplash -r  "addpath(genpath('matlab')), results = runtests('tests'), exit(0);"
   ```

In GNU Octave the same tests can be started with 
```
oruntests('tests')
```

Individual test scripts are also runnable as normal scripts, which 
makes them easier to debug in GNU Octave
