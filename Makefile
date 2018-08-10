external/Infomap: external/infomap/Infomap
	cp $< $@

external/infomap/Infomap:
	$(MAKE) -C external/infomap

# This target needs to be phony so it is run every time because only the other
# makefile can determine that there's nothing to be done.
.PHONY: external/infomap/Infomap
