default:
	@echo Hello.

compile:
	mkdir -p ./UV/
	g++ -std=c++11 -O4 -fopenmp ./rank.cpp -o ./UV/rank
	@echo 

clear:
	rm -f ./UV/rank

test:
	@echo Testing the rank routine.
	@echo Expect rank = 29876
	@echo Wait for it...
	@echo test_J2I_example.txt | time ./UV/rank 
