PROGRAM       = JToyFlowMain

version       = JTKT
CXX           = g++
CXXFLAGS      = -g -Wall   -Wno-deprecated -D$(version) 
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

CXXFLAGS  += $(INCS)

HDRSDICT = JToyFlowHistos.h 
           
HDRS	+= $(HDRSDICT)	JToyFlowInputs.h \
                        nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) JConst.h $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
