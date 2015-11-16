instdir = $(DESTDIR)$(docdir)
EXTRA_DIST = $(static_docs)

# Generate specification PDF documentation?
###########################################
if BUILD_LATEX

EXTRA_DIST += $(LATEX_BASENAME).pdf

LATEX_BASENAME = specification
LATEX_FILENAME = $(LATEX_BASENAME).tex

## COMMANDS:
PDFLATEX_CMD = pdflatex -interaction=nonstopmode
BIBTEX_CMD = bibtex


$(LATEX_BASENAME).pdf: $(LATEX_BASENAME).tex
	${PDFLATEX_CMD} $(LATEX_FILENAME)
	-${BIBTEX_CMD} $(LATEX_BASENAME)
	${PDFLATEX_CMD} $(LATEX_FILENAME)
	${PDFLATEX_CMD} $(LATEX_FILENAME)

endif

# Generate local API documentation?
###########################################
if BUILD_DOCS

DOXY_STAMP = doxygen-build.stamp
DOXY_FILE = doxy.conf
DOXY_AUX = header.html footer.html customstyles.css matlab.css DoxygenLayout.xml
HTML_AUX = $(top_srcdir)/doc/html/demoWriteGre.html $(top_srcdir)/doc/html/demoWriteRf.html $(top_srcdir)/doc/html/demoRead.html $(top_srcdir)/doc/html/Sequence.html
WEBSITE_AUX = $(LATEX_BASENAME).pdf logo_hires.png README.md favicon.ico images/gre.png images/rf2d.png images/overview.png

EXTRA_DIST += $(DOXY_FILE).in doxygen-build.stamp

#html/demoWrite.html: $(top_srcdir)/matlab/demoWrite.m
#	cd ../matlab && \
#	matlab -nodesktop -nosplash -nodisplay -r "make_matlab_doc('demoWrite.m','$(top_srcdir)/doc/html'); exit;" && \
#	cd ../doc

doxygen-build.stamp: $(DOXY_FILE) $(DOXY_AUX) $(HTML_AUX) $(top_srcdir)/src/*.h $(top_srcdir)/src/*.cpp *.dox
	doxygen $(DOXY_FILE)
	cp $(WEBSITE_AUX) html/
	touch doxygen-build.stamp

endif

doc_DATA = $(EXTRA_DIST)

install-data-local: $(EXTRA_DIST)
	$(mkinstalldirs) $(instdir)
	$(mkinstalldirs) $(instdir)/html
	if test -d html; then \
	  for f in html/*; do \
	    if test -f $$f; then \
		  $(INSTALL_DATA) $$f $(instdir)/html/; \
	    fi \
	  done; \
	  for dir in html/*; do \
	    if test -d $$dir; then \
	      b=`basename $$dir`; \
	      $(mkinstalldirs) $(instdir)/html/$$b; \
	      for f in $$dir/*; do \
	        if test -f $$f; then \
		      $(INSTALL_DATA) $$f $(instdir)/html/$$b; \
	        fi \
	      done \
	    fi \
	  done \
	fi; \
	if test -d $(LATEX_BASENAME).pdf; then \
	  $(INSTALL_DATA) $(LATEX_BASENAME).pdf $(instdir)/; \
	fi;

uninstall-local:
	rm -rf $(instdir)

clean-local:
	if test -d html; then rm -rf html; fi
	if test -f doxy-build.stamp; then rm -f doxy-build.stamp; fi
	if test -f doxy.conf; then rm -f doxy.conf; fi
	rm -f *.log *.aux *.synctex *.bbl *.blg *.toc *.out

install-html-local: $(EXTRA_DIST)
	cp -r html/* ../../pulseq_website/
