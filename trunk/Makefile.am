AUTOMAKE_OPTIONS = foreign 
ACLOCAL = aclocal
AUTOMAKE= automake
dist_doc_DATA = README INSTALL COPYING
clean-local:
	rm -rf autom4te.cache
SUBDIRS= \
	src bin

## Clean out SVN files and root dictionaries
dist-hook:
	rm -rf `find $(distdir) -name ".svn"` `find $(distdir) -name \*Dict\*`
