CC      = ifx
CFLAGS  = -O3 -march=native -ffast-math -fopenmp 

EXEC   = mol2_amarrage
OBJS   = lecture_mol2.o chargeur_covalence.o mol2_amarrage.o

all: $(EXEC)

run: $(EXEC)
	-@./$(EXEC) CoV_radii ligand.mol2 site.mol2

$(EXEC): $(OBJS)
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) $(CFLAGS) -o $@ $+

%.o: %.f90
	-@echo ""
	-@echo "Generating $@"
	-@$(CC) $(CFLAGS) -c $<


###------------------------------
### Cleaning
###------------------------------------------------------------

clean:
	-@rm -rf *.o
	-@rm -rf *.mod
	-@rm -rf $(EXEC)

clean_all: clean cleanSource
	-@rm -rf $(EXEC)

cleanSource:
	-@find . \( -name "*~" -o -name "*.old" -o -name "#*" \) -print -exec rm \{\} \;


.PHONY:  $(EXEC) clean clean_all cleanSource
