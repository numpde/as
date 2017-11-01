default:
	@echo Hello!

mark_unversioned:
	@echo -----------
	find -name 'UV' -type d | sed 's|.*|&/unversioned|' | xargs -L1 touch
	find -name 'unversioned' | xargs git add
	@echo -----------
	git status

fill_unversioned:
	for f in $$(find -name 'unversioned' -type f); do \
		d=$$(dirname $$f); \
		echo $$d; \
		ls $$d | grep -v '^unversioned$$' > $$f; \
		cat $$f; \
	done
	
	git commit -m "(ls unversioned directories)" $$(find -name 'unversioned' -type f)
