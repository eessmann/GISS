## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAVE_MODULES
%.c : %.mod
	@echo "/* $@" > $@
	@echo " * This is a generated file.  Please modify \`$<'" >> $@
	@echo " */" >> $@
	@echo "#line 1 \""$<"\"" >> $@
	@cat $< >> $@
else
%.c : %.mod
	@echo "/* $@" > $@
	@echo " * This is a generated file.  Please modify \`$<'" >> $@
	@echo " */" >> $@
	@echo "#line 1 \""$<"\"" >> $@
	@replace="s/g_module_check_init/gfs_init_`basename $< .mod`/g";\
	 sed $$replace < $< >> $@
endif

clean-generic:
	$(RM) $(BUILT_SOURCES)

dist-hook:
	for file in $(BUILT_SOURCES); do $(RM) $(distdir)/$$file; done
