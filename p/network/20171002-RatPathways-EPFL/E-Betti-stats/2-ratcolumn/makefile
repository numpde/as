
default:
	@echo Hi.

outdir:
	mkdir -p OUTPUT

compile:
	cd cpp; make compile

exp1: outdir compile
	python3 exp1*.py

exp2a: outdir
	F=$$(sed -n 's|.*\(OUTPUT.*-2a-.*.pkl\).*|\1|gp' exp2a*.py); \
	[ ! -e $$F ] && python3 exp2a*.py || python3 plot_exp2a*.py

exp2b: outdir
	F=$$(sed -n 's|.*\(OUTPUT.*-2b-.*.pkl\).*|\1|gp' exp2b*.py); \
	[ ! -e $$F ] && python3 exp2b*.py || python3 plot_exp2b*.py

