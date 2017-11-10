default:
	@echo Hello.

compile:
	mkdir -p ./UV/
	g++ -std=c++11 -O5 ./rank.cpp -o ./UV/rank

clear:
	rm -f ./UV/rank

test:
	@echo Expect rank = 
	@echo Wait for it...
	@echo test_J2I_example.txt | ./UV/rank 
