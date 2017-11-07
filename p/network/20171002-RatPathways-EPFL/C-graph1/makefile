default:
	@echo Hello!

step1:
	SCRIPT=step1_cliques.py; python3 $$SCRIPT  > $$SCRIPT.out

step2:
	SCRIPT=step2_laplacian.py; python3 $$SCRIPT  > $$SCRIPT.out

step3:
	SCRIPT=step3_1_eig.m; /NA/R2012a/bin/matlab -nosplash  < $$SCRIPT  > $$SCRIPT.out
