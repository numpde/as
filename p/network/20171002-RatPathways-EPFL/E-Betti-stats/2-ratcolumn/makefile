
default:
	@echo Hi.

outdir:
	mkdir -p OUTPUT

compile:
	cd cpp && make compile

exp1: outdir compile
	F=$$(ls exp1*.py); \
		( python3 $$F 1> OUTPUT/$$F.out.tmp 2> OUTPUT/$$F.err.tmp )
		&&
		( cd OUTPUT; mv $$F.out.tmp $$F.out.log; mv $$F.err.tmp $$F.err.log )

exp%: outdir
	F=$$(sed -n 's|.*\(OUTPUT.*$@.*.pkl\).*|\1|gp' $@*.py); \
		[ ! -e $$F ] && python3 $@*.py || python3 plot_$@*.py
