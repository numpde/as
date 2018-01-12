
default:
	@echo Hello.

0:
	mkdir -p OUTPUT/0_select/UV/
	python3 0_select.py 1> OUTPUT/0_select/info.log

p4:
	@echo Copying to public_html
	dest=~/public_html/files/20171212; \
	mkdir -p $$dest; \
	cp -v 4_*_a.html $$dest/index.html; \
	rsync -Rv $$(ls OUTPUT/4_*_a/*/*/*.png) $$dest;

6b:
	S=6_low-high-end_b.py; echo "export PYTHONUNBUFFERED=1; source ~/.bashrc; python3 $$S 1> $$S-cout.log 2> $$S-cerr.log" | at now
