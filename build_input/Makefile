NAME=myco

EVAL_CUT=0.001
NUM_THREADS=4
SEARCHSP=1000000000
OUTPUT_SIZE=3000
BLASTOPT_OLD=-Y $(SEARCHSP) -v $(OUTPUT_SIZE) -b $(OUTPUT_SIZE) -m 8 -e $(EVAL_CUT) -a $(NUM_THREADS)
BLASTOPT=-searchsp $(SEARCHSP) -max_target_seqs $(OUTPUT_SIZE) -evalue $(EVAL_CUT) -outfmt 6 -num_threads $(NUM_THREADS)

DOMCLUST_OPT=
OSNAME=$(shell uname -m| sed 's/  */_/')
BINDIR=../bin

## Comment out to use blast+; Uncomment out to use old blast
#OLDBLAST=true

all:	$(NAME).hom $(NAME).gene

prepare_blast:	$(NAME).pin

$(NAME).pin:	$(NAME)
ifeq ($(OLDBLAST),true)
	formatdb -i $<
else
	makeblastdb -in $< -dbtype prot
endif

.SECONDARY:	$(NAME).pin

$(NAME).blast:	$(NAME) $(NAME).pin
ifeq ($(OLDBLAST),true)
	blastall -p blastp -d $(NAME) -i $(NAME) $(BLASTOPT_OLD) -o $@
else
	blastp -db $(NAME) -query $(NAME) $(BLASTOPT) -out $@
endif

clustering:	$(NAME).clustout $(NAME).clusttree

$(NAME).hom:	$(NAME).blast
	./blast2homfile.pl -distconv $< > $@

$(NAME).gene:	$(NAME)
	./fasta2genefile.pl $<


$(NAME).clustout:	$(NAME).gene $(NAME).hom
	$(BINDIR)/domclust $(DOMCLUST_OPT) $(NAME).hom $(NAME).gene -ne1 | ../Script/addtit.pl -titfile=$(NAME).tit > $(NAME).clustout

$(NAME).clusttree:	$(NAME).gene $(NAME).hom
	$(BINDIR)/domclust $(DOMCLUST_OPT) $(NAME).hom $(NAME).gene -o1 -ne1 | ../Script/addtit.pl -titfile=$(NAME).tit -format=tree > $(NAME).clusttree

cleantmp:
	-rm *.blast *.pin *.psq *.phr *.log 2> /dev/null

clean:
	-rm *.gene *.hom *.blast *.pin *.psq *.phr *.tit *.log *.clustout *.clusttree 2> /dev/null
