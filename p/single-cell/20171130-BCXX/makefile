
default:
	@echo Hello.

0:
	mkdir -p OUTPUT/0_select/UV/
	python3 0_*.py 1> OUTPUT/0_select/info.log

p4:
	@echo Copying to public_html
	dest=~/public_html/files/20171212; \
	mkdir -p $$dest; \
	cp -v 4_*_a.html $$dest/index.html; \
	rsync -Rv $$(ls OUTPUT/4_*_a/*/*/*.png) $$dest;

