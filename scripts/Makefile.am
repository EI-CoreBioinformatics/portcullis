
install-exec-hook: 
	python3 setup.py install

uninstall-hook:
	pip3 uninstall junctools

junctools_scripts = \
	junctools/compare.py \
	junctools/convert.py \
	junctools/junction.py \
	junctools/markup.py \
	junctools/performance.py \
	junctools/set.py \
	junctools/split.py

dist-hook:
	cp $(junctools_scripts) $(distdir)/junctools