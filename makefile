#the compiler
CC = g++

#compiler flags 
# needs to be version 11 for boost to work
CFLAGS = -std=c++11
#include directory of some other codes
IDIR = ./commoncodes
PDIR = -I ../packages
INCDIR = -I ./commoncodes/ -I ./codefiles $(PDIR)
#codefile directory
CDIR = ./codefiles
FDIR = ./frontcodefiles
FTDIR = ${FDIR}/tests
#make directory
MKDIR_P = mkdir -p
ODIR = ./objects

TDIR = ./Tests
TDIRs = ${TDIR}/logs ${TDIR}/RunFiles/Ideal ${TDIR}/RunFiles/Ideal_Tab

#group codefiles
backfiles = $(addprefix $(CDIR)/, Jacob_Maker.hpp adiabatnab.hpp converge.hpp StructStruct.hpp PhysicalConsts.h EquationOfState.hpp )
#group objects
backobjects = $(addprefix $(ODIR)/, converge.o Jacob_Maker.o )
frontobjects = $(addprefix $(ODIR)/, TimeStructure.o GeneralisedStructure.o IrrEnvStructure2.o Struct_test_boost.o Ideal.o Tab_Ideal.o Poly_Time_test_boost.o Struct_test_boost_Irr.o Seager.o Stix.o Mix.o Abe.o StixStructure.o SeaStructure.o Anal_BC.o Tab_BC.o)

#executables
tests = Struct_test_boost
allexes = $(Idealexes) $(Idealexes_tab) $(Rockyexes)

all: directories test_directory $(Idealexes) $(Idealexes_tab) $(tests) $(Rockyexes)

test: test_directory $(tests) $(addprefix run-,$(tests))

directories: 
	$(MKDIR_P) ${ODIR}
test_directory:
	$(MKDIR_P) ${TDIRs}

Struct_test_boost: $(ODIR)/Struct_test_boost.o $(backobjects) $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o Struct_test_boost $(ODIR)/Struct_test_boost.o $(ODIR)/Ideal.o $(backobjects)

$(ODIR)/Ideal.o: $(CDIR)/IdealEquationOfState.cpp $(backfiles) $(CDIR)/EquationOfStateIdeal.hpp 
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Tab_Ideal.o: $(CDIR)/Tab_Ideal_EoS.cpp $(backfiles) $(IDIR)/2Dinterpolator.hpp $(IDIR)/File_error.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/converge.o: $(CDIR)/converge.cpp $(backfiles) $(CDIR)/block_tridiag_solve.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Jacob_Maker.o: $(CDIR)/Jacob_Maker.cpp $(backfiles)
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Struct_test_boost.o: $(FTDIR)/Struct_test_boost.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

run-%: %
	-./$^ > ./Tests/logs/$^.log

.PHONY : clean
clean:
	$(RM) $(backobjects) $(frontobjects) $(allexes) $(tests) $(addsuffix /* , ${TDIRs})

