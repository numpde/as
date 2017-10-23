default:
	@echo Hello.

connect:
	python3 extract.py > ratcolumn.txt
	grep    "^#" ratcolumn.txt > ratcolumn_head.txt
	grep -v "^#" ratcolumn.txt > ratcolumn_data.txt
	rm ratcolumn.txt

