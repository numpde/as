
default:
	@echo Type: make all

all:
	make clear
	make step1
	make step2
	make step3

clear:
	rm -f *.eps *.png *.pkl
	rm -rf __pycache__
	cd cpp; make clear

step1:
	@f=../D-cliques/cpp/variants/rank-9.cpp; [ -e $$f ] && cp -v $$f ./cpp/rank_localcopy.cpp
	@cd cpp; make compile

step2:
	@mkdir -p OUTPUT
	@python3 gnpbetti_compute.py

step3:
	@mkdir -p OUTPUT
	@python3 gnpbetti_plot.py
