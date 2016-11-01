AUTOMAKE_OPTIONS = subdir-objects

install-exec-hook: 
	python3 setup.py install

uninstall-hook:
	pip3 uninstall -y junctools

EXTRA_DIST = \
	junctools/compare.py \
	junctools/convert.py \
	junctools/junction.py \
	junctools/markup.py \
	junctools/performance.py \
	junctools/set.py \
	junctools/split.py

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf dist build junctools.egg-info