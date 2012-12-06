# sample Makefile.
# It compiles every .cpp files in the src/ directory to object files in obj/ directory, and build the ./run executable.
# It automatically handles include dependencies.

# You can modify it as you want (add libraries, compiler flags...).
# But if you want it to properly run on our benchmark service, don't rename or delete variables.

# using icc :
#COMPILER ?= $(ICC_PATH)icpc
# using gcc :
COMPILER ?= $(GCC_PATH)g++

# using icc :
#FLAGS ?= -std=c++0x -U__GXX_EXPERIMENTAL_COMPILER0X__ -xHOST -fast -w1 $(ICC_SUPPFLAGS)
# using gcc :
FLAGS ?= -std=c++0x -g -O3 -Wall $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
LDLIBS = 
#example if using Intel® Threading Building Blocks :
#LDLIBS = -ltbb -ltbbmalloc

EXECUTABLE = run

TEAM_ID = 30f641898926b29c74b73e3275bee3e3 # put your 32chars team id here and you will be able to submit your program from command line using "make submit" 

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:src/%.cpp=obj/%.o)

all: release

release: objdir $(OBJS)
	$(COMPILER) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

objdir:
	mkdir -p obj

obj/%.o: src/%.cpp
	$(COMPILER) $(FLAGS) -o $@ -c $<


zip:
ifdef TEAM_ID
	zip $(strip $(TEAM_ID)).zip -r Makefile src
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif

submit: zip
ifdef TEAM_ID
	curl -F "file=@$(strip $(TEAM_ID)).zip" -L http://www.intel-software-academic-program.com/contests/ayc/2012-11/upload/upload.php
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif

clean:
	rm -f obj/* run
dist-clean: clean
	rm -f $(EXECUTABLE) *~ .depend *.zip
	
#automatically handle include dependencies
depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	@$(foreach SRC, $(SRCS), $(COMPILER) $(FLAGS) -MT $(SRC:src/%.cpp=obj/%.o) -MM $(SRC) >> .depend;)

run-tests: release
	rm -f tests/*/*.out
	for TEST in $(TESTS); do \
		tput setaf 4; \
		echo Running $$TEST; \
		tput sgr 0; \
		echo /usr/bin/time --output=tests/$$TEST/time.out -f "%E real, %U user, %S sys" \
			tests/$$TEST/script.sh \
			-flights tests/$$TEST/flights.txt \
			-alliances tests/$$TEST/alliances.txt \
			-work_hard_file tests/$$TEST/work_hard.out \
			-play_hard_file tests/$$TEST/play_hard.out ; \
		/usr/bin/time --output=tests/$$TEST/time.out -f "%E real, %U user, %S sys" \
			tests/$$TEST/script.sh \
			-flights tests/$$TEST/flights.txt \
			-alliances tests/$$TEST/alliances.txt \
			-work_hard_file tests/$$TEST/work_hard.out \
			-play_hard_file tests/$$TEST/play_hard.out ; \
		tput setaf 2; \
		cat tests/$$TEST/time.out; \
		tput setaf 1; \
		diff -u tests/$$TEST/work_hard.txt tests/$$TEST/work_hard.out; \
		diff -u tests/$$TEST/play_hard.txt tests/$$TEST/play_hard.out; \
		tput sgr0; \
	done

unit-tests:
	mkdir -p obj
	g++ $(FLAGS) src/test_UniqueId.cc -o obj/test_UniqueId && obj/test_UniqueId
	g++ $(FLAGS) src/test_Alliances.cc -o obj/test_Alliances && obj/test_Alliances

test:
	make run-tests TESTS="correctness1 scalability1 scalability2 scalability3"

include .depend	
