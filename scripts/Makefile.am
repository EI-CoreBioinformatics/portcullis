AUTOMAKE_OPTIONS = subdir-objects

install-exec-hook:
	export PYTHONPATH=$(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	mkdir -p $(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	cd junctools && \
	python3 setup.py install --prefix=$(prefix) --record install_files.txt && \
	cd ../portcullis && \
	python3 setup.py install --prefix=$(prefix) --record install_files.txt && \
	cd ..


uninstall-hook:
	cat junctools/install_files.txt | xargs rm -rf	
	cat portcullis/install_files.txt | xargs rm -rf	

EXTRA_DIST = \
	junctools/junctools/compare.py \
	junctools/junctools/convert.py \
	junctools/junctools/junction.py \
	junctools/junctools/markup.py \
	junctools/junctools/performance.py \
	junctools/junctools/set.py \
	junctools/junctools/split.py \
	junctools/junctools/gtf.py \
	junctools/junctools/__init__.py \
	junctools/junctools/__main__.py \
	portcullis/portcullis/rule_filter.py \
	portcullis/portcullis/__init__.py

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf junctools/dist junctools/build junctools/junctools.egg-info
	-rm -rf portcullis/dist portcullis/build portcullis/junctools.egg-info
