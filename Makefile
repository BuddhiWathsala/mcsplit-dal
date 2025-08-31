all:
	g++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -pthread -o bin/run.o
build:
	g++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -pthread -o bin/run.o
linuxbuild:
	g++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -g -pthread -o bin/run.o
macbuild:
	/usr/bin/clang++ -O3 -ffast-math -pedantic -Wall -std=c++17 *.cpp -pthread -g -o bin/run.o -I/opt/homebrew/Cellar/boost/1.86.0_2/include/ -I/opt/homebrew/Cellar/argp-standalone/1.3/include -L/opt/homebrew/Cellar/argp-standalone/1.3/lib/ -largp
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
	@echo 'Solution (case 16) = 51'
	./bin/run.o min_max ./data/tests/case17/pattern ./data/tests/case17/target -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 17) = 40'
	./bin/run.o min_max ./data/tests/case18/pattern ./data/tests/case18/target -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 18) = 40'
	./bin/run.o min_max ./data/tests/case19/g22 ./data/tests/case19/g57 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 19) = 8'
	./bin/run.o min_max ./data/tests/case20/g5 ./data/tests/case20/g61 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 20) = 10'
	./bin/run.o min_max ./data/tests/case21/g9 ./data/tests/case21/g80 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 21) = 12'
	./bin/run.o min_max ./data/tests/case22/pattern ./data/tests/case22/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case23/g8 ./data/tests/case23/g37 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 23) = 30'
	./bin/run.o min_max ./data/tests/case24/pattern1 ./data/tests/case24/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case25/pattern ./data/tests/case25/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case26/g27 ./data/tests/case26/g5 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 26) = 11'
	./bin/run.o min_max ./data/tests/case27/g44 ./data/tests/case27/g6 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 27) = 19'
	./bin/run.o min_max ./data/tests/case28/g5 ./data/tests/case28/g27 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 28) = 11'
	./bin/run.o min_max ./data/tests/case29/g6 ./data/tests/case29/g44 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 29) = 19'
	./bin/run.o min_max ./data/tests/case30/g13 ./data/tests/case30/g6 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 30) = 15'
	./bin/run.o min_max ./data/tests/case31/g25 ./data/tests/case31/g6 -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 31) = 13'
	./bin/run.o min_max ./data/tests/case32/pattern ./data/tests/case32/target -l -c -q -t ${MAX_TIME}
	@echo 'Solution (case 32) = 120'
	./bin/run.o min_max ./data/tests/case33/pattern ./data/tests/case33/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case34/g3 ./data/tests/case34/g50 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case35/pattern ./data/tests/case35/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case36/g4 ./data/tests/case36/g10 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case37/g7 ./data/tests/case37/g13 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case38/pattern ./data/tests/case38/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case39/pattern ./data/tests/case39/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case40/pattern ./data/tests/case40/target -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case41/g15 ./data/tests/case41/g57 -l -c -q -t ${MAX_TIME}
	./bin/run.o min_max ./data/tests/case42/g6 ./data/tests/case42/g17 -l -c -q -t ${MAX_TIME}
