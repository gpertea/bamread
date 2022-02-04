GDIR := ../gclib

## assumed htslib has been pulled from https://github.com/gpertea/htslib
HTSLIB_SRC := ../htslib
HTSLIB := $(abspath $(HTSLIB_SRC))

ifeq ($(wildcard $(HTSLIB)),)
  $(info Clone https://github.com/gpertea/htslib in the parent directory)
  $(error $(HTSLIB_SRC) source missing)
#else
#  $(info $(HTSLIB_SRC) found!)
endif

#my branch of htslib includes libdeflate:
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a
LIBBZ2 := ${HTSLIB}/xlibs/lib/libbz2.a
INCDIRS := -I. -I${GDIR} -I${HTSLIB}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    



CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS  := -Wall -Wextra ${INCDIRS} $(MARCH) \
 -D_REENTRANT -std=c++11 -fno-strict-aliasing -fno-exceptions -fno-rtti

ifneq (,$(filter %release %static %static-cpp, $(MAKECMDGOALS)))
  #release build
  CFLAGS := -g -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
else
  CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${HTSLIB}
endif

ifneq ($(findstring static,$(MAKECMDGOALS)),) 
 # static or static-cpp found
 ifneq ($(findstring static-cpp,$(MAKECMDGOALS)),) 
    #not a full static build, only c/c++ libs
    LDFLAGS := -static-libgcc -static-libstdc++ ${LDFLAGS}
 else
    #full static build
    LDFLAGS := -static -static-libgcc -static-libstdc++ ${LDFLAGS}
 endif
endif

%.o : %.cpp
	${CXX} ${CFLAGS} -c $< -o $@

# C/C++ linker
LINKER  := $(if $(LINKER),$(LINKER),g++)

LIBS := ${HTSLIB}/libhts.a ${LIBLZMA} ${LIBDEFLATE} ${LIBBZ2} -lz -lm -lpthread
OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o \
        GSam.o

# Compiling for Windows with MinGW/MSYS2?
ifneq ($(findstring -mingw,$(shell $(CXX) -dumpmachine 2>/dev/null)),)
 LIBS += -lregex -lws2_32
endif

.PHONY : all
all release static static-cpp debug: samread

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


