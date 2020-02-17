GCLDIR := ../gclib
HTSLIB := ../htslib
LIBDEFLATE := /ccb/sw/lib/libdeflate.a

SEARCHDIRS := -I. -I${GCLDIR} -I${HTSLIB}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

CC      := g++

BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) \
 -D_REENTRANT -fno-strict-aliasing -fno-exceptions -fno-rtti

#add the link-time optimization flag if gcc version > 4.5


ifeq ($(findstring release,$(MAKECMDGOALS)),)
  CFLAGS := -g -DDEBUG -DGDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
else
  CFLAGS := -g -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := ${HTSLIB}/libhts.a ${LIBDEFLATE} -llzma -lbz2 -lz -lm -lcurl -lcrypto -lpthread
OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GStr.o \
        GSam.o
 
.PHONY : all
all:    bamread
debug:  bamread

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/GBase.h 
bamread.o : ./GSam.h
GSam.o : ./GSam.h

bamread: $(OBJS) bamread.o
	${LINKER} ${LDFLAGS} $(GCC45OPTS) $(GCC45OPTMAIN) -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} bamread bamread.o* bamread.exe $(OBJS)
	@${RM} core.*


