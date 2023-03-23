SYSNAME := $(shell uname -m | sed 's/  */_/')
#BINDIR := bin-$(SYSNAME)
#BINDIR2 := bin-$(SYSNAME).2
#BINDIR3 := bin-$(SYSNAME).3

BINDIR := bin
BINDIR2 := bin.2
BINDIR3 := bin.3

#sinclude ../conf/$(SYSNAME).conf

CC=gcc
#CC=icc
CFLAGS0=-w -fopenmp
CFLAGS=-O3 $(CFLAGS0)
#CFLAGS += -DEXPERIMENTAL

#CFLAGS=-O2 -g $(CFLAGS0)
#CFLAGS=-g $(CFLAGS0)

TESTOPT=-S60 -C80 -V.6 -ao0.8
#TESTOPT=-S60 -C80 -V.6 -ao0.8 -Cl.2
#TESTOPT=-S60 -C80 -V.6 -ao0.8 -Omincutcnt=2

OBJS := domclust.o util.o graph.o memalloc.o namehash.o hash.o \
	cluster.o clustout.o collectclust.o vararray.o plist.o bin.o seqreg.o \
	spec.o sptree.o domain.o readfile.o dump.o neighbor.o \
	postproc.o mergeclust.o newali.o djset.o clusttree.o taxmap.o sstack.o

OBJS2 = $(OBJS) genedp.o domcut.o calnewali.o mergehit.o poslist.o hitlist.o

SELOBJS = select.o memalloc.o plist.o hash.o namehash.o vararray.o
SPINDEXOBJS = spindex.o
SLINKOBJS = slink.o memalloc.o plist.o hash.o namehash.o

all2 : CFLAGS += -DWITH_NEIGHBOR -DWITH_DOMCUT -DWITH_CALNEWALI

all_targetN : CFLAGS += -DLARGE

all:	$(BINDIR)
	cd $(BINDIR); \
	$(MAKE) -f ../Makefile VPATH=.. SYSNAME=$(SYSNAME) all_target


all2:	$(BINDIR2)
	cd $(BINDIR2); \
	$(MAKE) -f ../Makefile VPATH=.. SYSNAME=$(SYSNAME) all_target2

all3:	$(BINDIR3)
	cd $(BINDIR3); \
	$(MAKE) -f ../Makefile VPATH=.. SYSNAME=$(SYSNAME) all_targetN

all_target:	domclust select spindex blastconv

all_target2:	domclust2

all_targetN:	domclustN

$(BINDIR):
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi

$(BINDIR2):
	@if [ ! -d $(BINDIR2) ]; then mkdir $(BINDIR2); fi

$(BINDIR3):
	@if [ ! -d $(BINDIR3) ]; then mkdir $(BINDIR3); fi


domclust:	$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

domclust2:	$(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) -o $@ -lm

domclustN:	$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

sptree:	sptree.c plist.c namehash.c memalloc.c hash.c spec.c
	$(CC) -g sptree.c plist.c namehash.c memalloc.c hash.c spec.c -DDEBUGMAIN_SPTREE -o $@

clean:
	-rm -rf $(BINDIR)
clean2:
	-rm -rf $(BINDIR2)

select:	$(SELOBJS)
	$(CC) $(CFLAGS) $(SELOBJS) -lgdbm -o $@

spindex:	$(SPINDEXOBJS)
	$(CC) $(CFLAGS) $(SPINDEXOBJS) -lgdbm -o $@

blastconv:	blastconv.o
	$(CC) $(CFLAGS) blastconv.o -o $@

test:	all
	cd tst;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test.in | \
		../Script/cmpr.pl -OVLPRATIO2=.2 -diff test.out -;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test2.in | \
		../Script/cmpr.pl -OVLPRATIO2=.2 -diff test2.out -;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test3.in | \
		../Script/cmpr.pl -OVLPRATIO2=.2 -diff test3.out -

testout:	all
	cd tst;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test.in > test.out;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test2.in > test2.out;\
	../$(BINDIR)/domclust stdin -S $(TESTOPT) < test3.in > test3.out

dist:	distrib/domclust
	for f in $(OBJS:.o=.c) $(OBJS:.o=.h); do \
		if [ -f $$f ]; then \
			Script/distfilter.pl $$f > distrib/domclust/$$f; \
		fi \
	done
	cd distrib; Proc.build

distrib/domclust:
	-mkdir -p distrib/domclust

tarfile:
	tar cvf domclust.tar *.c *.h Makefile Script conf tst
