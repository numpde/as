default:
	@echo Hello.
	@echo Try: make V=3 compile test

compile:
	g++ -std=c++11 -O5 ./rank-${V}.cpp -o ./rank
	@echo 

clean:
	rm -f ./rank

test:
	@echo Testing the rank routine.
	@echo Expect rank = 29876
	@echo Wait for it...
	@echo ../test_J2I_example.txt | time ./rank 
