SHELL = /bin/sh
TARGET = libgtmcmc_wb.a
CSRCS = tmcmc_wb.cpp
DEPCSRC = $(CSRCS)
OBJS = $(CSRCS:.cpp=.o)
INCFILES = tmcmc_wb.h
INCDIR=../include
INCL = -I. 
CPP = mpic++
CPPFLAGS = -O2 -std=gnu++11 -ggdb $(INCL)
RANLIB = ranlib
AR = ar cr

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)
	for i in $(INCFILES) ;\
	  do /bin/rm -f $(INCDIR)/$${i} ;\
	  /bin/cp $$i $(INCDIR)/$${i} ;\
	done

clean:
	rm -f $(OBJS) $(TARGET)

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $*.cpp
