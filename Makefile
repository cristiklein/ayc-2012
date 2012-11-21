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
FLAGS ?= -std=c++0x -O3 -Wall $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
LDLIBS = 
#example if using Intel� Threading Building Blocks :
#LDLIBS = -ltbb -ltbbmalloc

EXECUTABLE = run

TEAM_ID = 30f641898926b29c74b73e3275bee3e3 # put your 32chars team id here and you will be able to submit your program from command line using "make submit" 

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:src/%.cpp=obj/%.o)

all: release

release: $(OBJS)
	$(COMPILER) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

obj/%.o: src/%.cpp
	$(COMPILER) $(FLAGS) -o $@ -c $<


zip: dist-clean
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

test: release
	rm -f scenarios/newinput/play_hard.txt
	rm -f scenarios/newinput/work_hard.txt
	scenarios/newinput/script.sh
	diff -u scenarios/newinput/play_hard.txt.expected scenarios/newinput/play_hard.txt
	diff -u scenarios/newinput/work_hard.txt.expected scenarios/newinput/work_hard.txt
	rm -f scenarios/newinput/play_hard.txt
	rm -f scenarios/newinput/work_hard.txt

test2: release
	rm -f scenarios/newinput2/play_hard.txt
	rm -f scenarios/newinput2/work_hard.txt
	scenarios/newinput2/script.sh
	diff -u scenarios/newinput2/play_hard.txt.expected scenarios/newinput2/play_hard.txt
	diff -u scenarios/newinput2/work_hard.txt.expected scenarios/newinput2/work_hard.txt
	rm -f scenarios/newinput2/play_hard.txt
	rm -f scenarios/newinput2/work_hard.txt

include .depend	
