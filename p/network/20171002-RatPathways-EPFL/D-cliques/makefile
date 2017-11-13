
default:
	@echo Hi.

outdir:
	mkdir -p OUTPUT

compile:
	cd cpp; make compile

step1: outdir compile
	python3 step1*.py

exp2a: outdir
	python3 exp2a*.py

exp2b: outdir
	python3 exp2b*.py

