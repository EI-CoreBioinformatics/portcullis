AUTOMAKE_OPTIONS = subdir-objects

if MAKE_PY_NOINSTALL

install-exec-hook:
	mkdir -p $(DESTDIR)$(datadir)/portcullis/scripts/portcullis
	cp -r portcullis/portcullis/*.py $(DESTDIR)$(datadir)/portcullis/scripts/portcullis 


uninstall-hook:
	rm -rf $(DESTDIR)$(datadir)/portcullis/scripts

endif

if MAKE_PY_INSTALL

install-exec-hook:
	mkdir -p $(DESTDIR)$(datadir)/portcullis/scripts/portcullis
	cp -r portcullis/portcullis/*.py $(DESTDIR)$(datadir)/portcullis/scripts/portcullis
	cd junctools && \
	$(PYTHON_INT_PATH) setup.py install --record install_files.txt && \
	cd ../portcullis && \
	$(PYTHON_INT_PATH) setup.py install --record install_files.txt && \
	cd ..


uninstall-hook:
	rm -rf $(DESTDIR)$(datadir)/portcullis/scripts
	cd junctools && \
	cat install_files.txt | xargs rm -rf && \
	cd ../portcullis && \
	cat install_files.txt | xargs rm -rf && \
	cd ..

endif

if MAKE_PY_INSTALL_CUSTOM


install-exec-hook:
	mkdir -p $(DESTDIR)$(datadir)/portcullis/scripts/portcullis
	cp -r portcullis/portcullis/*.py $(DESTDIR)$(datadir)/portcullis/scripts/portcullis
	mkdir -p $(DESTDIR)$(pythondir)
	export PYTHONPATH=$(DESTDIR)$(pythondir) && \
	cd junctools && \
	$(PYTHON_INT_PATH) setup.py install --prefix=$(DESTDIR)$(prefix) --record install_files.txt && \
	cd ../portcullis && \
	$(PYTHON_INT_PATH) setup.py install --prefix=$(DESTDIR)$(prefix) --record install_files.txt && \
	cd ..

uninstall-hook:
	rm -rf $(DESTDIR)$(datadir)/portcullis/scripts
	cd junctools && \
	cat install_files.txt | xargs rm -rf && \
	cd ../portcullis && \
	cat install_files.txt | xargs rm -rf && \
	cd ..


endif



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
