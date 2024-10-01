<div align="center">
<a href="http://pulseq.github.io/" target="_blank">
<img src="http://pulseq.github.io/logo_hires.png" width="250" alt="Pulseq logo"></img>
</a>
</div>

# Website

This repository stores the public website automatically generated from the *Pulseq* source code using doxygen and matlab. The resulting website is accessible here:

http://pulseq.github.io/

The source code is here: https://github.com/pulseq/pulseq/doc

the workflow: 
1: at the parent level run 
    make install-html
  this will create a directory named pulseq_website on the lvel above the partn and will copy all the necessary files over

2: check the webpage locally 

3: copy over to the checked-out location of http://pulseq.github.io/ and commit

Optional step 0: prior to the above you may wan to run: 
    matlab ./matlab/make_all_docs.m


