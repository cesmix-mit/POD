CPP=/data/scratch/teoc/projects/hlc/Exasim/External/llvm/prefix/bin/clang++
THB=/data/scratch/teoc/projects/hlc/Exasim/External/halide/build
THS=/data/scratch/teoc/projects/hlc/Exasim/External/halide/Halide
CILK=/data/scratch/teoc/projects/hlc/Exasim/External/llvm/prefix/lib/clang/12.0.0/lib/x86_64-unknown-linux-gnu/

hbin/gen: gen.cpp
	clang++ gen.cpp ${THS}/tools/GenGen.cpp -g -std=c++17 -fno-rtti -I${THB}/include -L${THB}/src -lHalide -lpthread -ldl -o hbin/gen -L${THB}/src -L${CILK}

pod1.o: hbin/gen
	LD_LIBRARY_PATH=${THB}/src:${CILK} hbin/gen -g pod1 -o . target=host -e static_library,stmt,stmt_html,c_header,object

snapshot.o: hbin/gen
	LD_LIBRARY_PATH=${THB}/src:${CILK} hbin/gen -g snapshot -o . target=host -e static_library,stmt,stmt_html,c_header,object

pod: main.cpp pod.cpp snapshot.o readinputfiles.cpp
	${CPP} -std=c++17 -O3 -lpng -ljpeg -lz -lpthread -ldl -lopencilk -lopencilk-personality-cpp -lopencilk-personality-c -I${THS}/tools -I${THB}/include -L/data/scratch/teoc/projects/hlc/Exasim/External/llvm/prefix/lib/clang/12.0.0/lib/x86_64-unknown-linux-gnu/ -Wall -llapack -lblas -g  main.cpp snapshot.o -o pod

default: main.cpp pod.cpp
	${CPP} -std=c++17 -O3 -Wall -llapack -lblas -g  main.cpp -o pod

all: pod
clean:
	rm pod hbin/gen pod1.h *.stmt *.html *.o *.a


run: pod
	cd /data/scratch/teoc/projects/pod/POD/fitting/Ta;\
	LD_LIBRARY_PATH=${THB}/src:${CILK} ./../../src/pod pod.txt data.txt
