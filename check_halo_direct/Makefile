CPP := mpicxx
CPPFLAGS := -O0
SRCDIR := ./src
OBJDIR := ./obj

HACC_SRC_PATH := /projects/exasky/HACC/trunk
HACC_EXE_PATH := /projects/exasky/HACC/trunk/darwin.cpu

INCLUDES := -I $(SRCDIR) -I$(HACC_SRC_PATH)/genericio -I$(HACC_SRC_PATH)/cosmotools/algorithms/halofinder -I$(HACC_SRC_PATH)/cosmotools/algorithms/mergertree_ct -I$(HACC_SRC_PATH)/cosmotools/common -I./src

#
OBJECTS := $(OBJDIR)/main.o $(OBJDIR)/particle_check.o  $(OBJDIR)/compare_dist.o $(OBJDIR)/match_pos.o $(OBJDIR)/sod_check.o $(OBJDIR)/halo_check.o $(OBJDIR)/mhn.o $(OBJDIR)/Halos_test.o


#linking
run_test: $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L$(HACC_EXE_PATH)/mpi/lib -lCosmoTools -lGenericIOMPI -o run_test -fopenmp


#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	mkdir -p obj
	mkdir -p outputs
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp
$(OBJDIR)/sod_check.o: $(SRCDIR)/sod_check.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/sod_check.cpp -o $(OBJDIR)/sod_check.o -fopenmp
$(OBJDIR)/halo_check.o: $(SRCDIR)/halo_check.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/halo_check.cpp -o $(OBJDIR)/halo_check.o -fopenmp
$(OBJDIR)/particle_check.o: $(SRCDIR)/particle_check.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/particle_check.cpp -o $(OBJDIR)/particle_check.o -fopenmp
$(OBJDIR)/Halos_test.o: $(SRCDIR)/Halos_test.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/Halos_test.cxx -o $(OBJDIR)/Halos_test.o -fopenmp
$(OBJDIR)/compare_dist.o: $(SRCDIR)/compare_dist.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/compare_dist.cpp -o $(OBJDIR)/compare_dist.o -fopenmp
$(OBJDIR)/match_pos.o: $(SRCDIR)/nonid_match.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/nonid_match.cpp -o $(OBJDIR)/match_pos.o -fopenmp
$(OBJDIR)/mhn.o: $(SRCDIR)/MurmurHashNeutral2.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -g -c $(SRCDIR)/MurmurHashNeutral2.cpp -o $(OBJDIR)/mhn.o -fopenmp




clean:
	rm $(OBJDIR)/*.o run_test
