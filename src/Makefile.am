bin_PROGRAMS = parsemr

if BUILD_TESTS
  TESTS = testparser.py
endif

SOURCES = ExternalSequence.cpp ExternalSequence.h
parsemr_SOURCES = $(SOURCES) parsemr.cpp 

EXTRA_DIST = $(TESTS)


clean-local:
	rm -f *.o 