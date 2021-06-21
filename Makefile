GDIR := ../gclib

## assumed htslib has been pulled from https://github.com/gpertea/htslib
HTSLIB := ../htslib
#my branch of htslib includes libdeflate:
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a
INCDIRS := -I. -I${GDIR} -I${HTSLIB}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

CC      := g++

BASEFLAGS  := -Wall -Wextra ${INCDIRS} $(MARCH) \
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
OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o \
        GSam.o

# Compiling for Windows with MinGW/MSYS2?
ifneq ($(findstring -mingw,$(shell $(CC) -dumpmachine 2>/dev/null)),)
 LIBS += -lregex -lws2_32
endif

.PHONY : all
all release static debug: samread

$(OBJS) : $(GDIR)/GBase.h $(GDIR)/GBase.h 
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


