
if BUILD_TESTS
if BUILD_MATLAB_TESTS
  TESTS = testmatlab.py
endif
endif

# Generate local API documentation?
if BUILD_DOCS


HTMLDIR = $(top_srcdir)/doc/html

doc_DATA = matlab_build.stamp

#SUFFIXES = .html .m
#$(HTMLDIR)/%.html: %.m

#.m.html:
#	matlab -nodesktop -nosplash -nodisplay -r "make_matlab_doc('$<','$(HTMLDIR)'); exit;"

matlab_build.stamp: $(HTMLDIR)/demoWriteGre.html $(HTMLDIR)/demoWriteRf.html $(HTMLDIR)/demoRead.html $(HTMLDIR)/Sequence.html
	touch matlab_build.stamp

$(HTMLDIR)/demoWriteGre.html: demoWriteGre.m make_matlab_doc.m
	matlab -wait -nodesktop -nosplash  -r "make_matlab_doc('$<','$(HTMLDIR)'); exit;"

$(HTMLDIR)/demoWriteRf.html: demoWriteRf.m make_matlab_doc.m
	matlab -wait -nodesktop -nosplash  -r "make_matlab_doc('$<','$(HTMLDIR)'); exit;"

$(HTMLDIR)/demoRead.html: demoRead.m make_matlab_doc.m
	matlab -wait -nodesktop -nosplash  -r "make_matlab_doc('$<','$(HTMLDIR)'); exit;"

$(HTMLDIR)/Sequence.html: +mr/@Sequence/Sequence.m make_matlab_doc.m
	matlab -wait -nodesktop -nosplash  -r "make_matlab_doc('$<','$(HTMLDIR)'); exit;"
endif


clean-local:
	if test -d html; then rm -rf html; fi
