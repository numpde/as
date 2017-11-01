default:
	@echo Hello!

step1:
	SCRIPT=step1-cliques.py; python3 $$SCRIPT  > $$SCRIPT.out

step2:
	SCRIPT=step2-laplacian.py; python3 $$SCRIPT  > $$SCRIPT.out

step3:
	SCRIPT=step3-1-eig.m; /NA/R2012a/bin/matlab -nosplash  < $$SCRIPT  > $$SCRIPT.out
