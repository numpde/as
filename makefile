default:
	@echo Hello!

mark_unversioned:
	@echo -----------
	find -name 'UV' -type d | sed 's|.*|&/unversioned|' | xargs -L1 touch
	find -name 'unversioned' | xargs git add
	@echo -----------
	git status
