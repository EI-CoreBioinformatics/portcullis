AUTOMAKE_OPTIONS = subdir-objects

install-exec-hook:
	PYTHONPATH=$(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	mkdir -p $(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	python3 setup.py install --prefix=$(prefix) --record install_files.txt

uninstall-hook:
	cat install_files.txt | xargs rm -rf	

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