all:
	g++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -pthread -o bin/run.o
build:
	g++ -O3 -ffast-math -g -pedantic -Wall -std=c++17 *.cpp -pthread -o bin/run.o
linuxbuild:
	g++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -g -pthread -o bin/run.o
macbuild:
	/usr/bin/clang++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -pthread -g -o bin/run.o -I/opt/homebrew/Cellar/argp-standalone/1.3/include -L/opt/homebrew/Cellar/argp-standalone/1.3/lib/ -largp
run:
	./bin/run.o min_max ./data/tests/case1/pattern ./data/tests/case1/target -l -c

MAX_TIME := 10
test:
	./bin/run.o min_max ./data/tests/case1/pattern ./data/tests/case1/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case2/pattern ./data/tests/case2/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case3/g26 ./data/tests/case3/g92 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case4/g17 ./data/tests/case4/g80 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case5/g27 ./data/tests/case5/g52 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case6/g50 ./data/tests/case6/g92 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case7/g11 ./data/tests/case7/g93 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case8/g10 ./data/tests/case8/g51 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case9/pattern ./data/tests/case9/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case10/pattern ./data/tests/case10/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case12/g10 ./data/tests/case12/g105 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case13/g10 ./data/tests/case13/g104 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case15/pattern ./data/tests/case15/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case16/pattern ./data/tests/case16/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case17/pattern ./data/tests/case17/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case18/pattern ./data/tests/case18/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case19/g22 ./data/tests/case19/g57 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case20/g5 ./data/tests/case20/g61 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case21/g9 ./data/tests/case21/g80 -l -c -q -t ${MAX_TIME}
