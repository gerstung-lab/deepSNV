all: roxygen build

roxygen: deepSNV
	#R CMD roxygen -u deepSNV
	find deepSNV/man  | xargs rm -rf
	R -e 'library(roxygen2); roxygenize("deepSNV")'
	#./patch.sh
	patch --no-backup-if-mismatch deepSNV/NAMESPACE patches/NAMESPACE.patch
	autoconf -o deepSNV/configure deepSNV/configure.ac
	
check: roxygen clean
	R CMD check deepSNV
	
build: roxygen clean
	R CMD build --no-build-vignettes deepSNV
	mv *.tar.gz builds
	
install: deepSNV
	R CMD install deepSNV
	
clean: deepSNV*
	find deepSNV -name '.*' | xargs rm -rf
	find deepSNV -name '*.o' | xargs rm -rf
	find deepSNV -name '*.so' | xargs rm -rf
	find deepSNV -name 'config.*' | xargs rm -rf
	find deepSNV -name 'symbols.rds' | xargs rm -rf
	find . -name 'deepSNV.Rcheck' | xargs rm -rf
	find . -name 'autom4te.cache' | xargs rm -rf
	
bioc: roxygen clean
	#find  bioC/devel/deepSNV -type f | grep -v '/\.' | xargs rm 
	#rsync -avp --exclude='.*' deepSNV/ bioc/devel/deepSNV/
	
