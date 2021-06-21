GCLDIR := ../gclib
## assumed htslib has been pulled from https://github.com/gpertea/htslib 
## (branch for_gclib) and build_lib.sh script was run in that directory
HTSLIB := ../htslib
#my branch of htslib includes libdeflate:
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a
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
 -D_REENTRANT -std=c++11 -fno-strict-aliasing -fno-exceptions -fno-rtti

#add the link-time optimization flag if gcc version > 4.5


ifeq ($(findstring release,$(MAKECMDGOALS)),)
  CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
else
  CFLAGS := -g -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := ${HTSLIB}/libhts.a ${LIBLZMA} ${LIBDEFLATE} -lbz2 -lz -lm -lpthread
OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GStr.o \
        GSam.o

# Compiling for Windows with MinGW?
ifneq ($(findstring -mingw,$(shell $(CC) -dumpmachine 2>/dev/null)),)
 LIBS += -lregex -lws2_32
endif

.PHONY : all
all release static debug: samread

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/GBase.h 
samread.o : ./GSam.h
GSam.o : ./GSam.h

${HTSLIB}/libhts.a: 
	cd ${HTSLIB} && ./build_lib.sh

samread: ${HTSLIB}/libhts.a $(OBJS) samread.o
	${LINKER} ${LDFLAGS} $(GCC45OPTS) $(GCC45OPTMAIN) -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} samread samread.o* samread.exe $(OBJS)
	@${RM} core.*


