GCLDIR := ../gclib
BAM_INCDIR := ./samtools-0.1.18
BAM_LIBDIR := ./samtools-0.1.18

SEARCHDIRS := -I. -I${GCLDIR} -I${BAM_INCDIR}

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
  CFLAGS := -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS := -L${BAM_LIBDIR}
else
  CFLAGS := -g -DDEBUG -DGDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${BAM_LIBDIR}
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := -lz -lbam
OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GStr.o \
        GXam.o
 
.PHONY : all
all:    bamread
debug:  bamread

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/GBase.h 
bamread.o : ./GXam.h
GXam.o : ./GXam.h

bamread: $(OBJS) bamread.o
	${LINKER} ${LDFLAGS} $(GCC45OPTS) $(GCC45OPTMAIN) -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} bamread bamread.o* bamread.exe $(OBJS)
	@${RM} core.*


