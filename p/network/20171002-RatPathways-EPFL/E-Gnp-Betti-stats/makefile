
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
	cd cpp; make compile

step2:
	python3 gnpbetti_compute.py

step3:
	python3 gnpbetti_plot.py
