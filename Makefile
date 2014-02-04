all: roxygen build

roxygen: deepSNV
	#R CMD roxygen -u deepSNV
	R -e 'library(roxygen2); roxygenize("deepSNV", unlink.target=TRUE)'
	#./patch.sh
	patch --no-backup-if-mismatch deepSNV/NAMESPACE patches/NAMESPACE.patch
	autoconf -o deepSNV/configure deepSNV/configure.ac
	
check: roxygen
	R CMD check deepSNV
	
build: roxygen
	R CMD build --no-build-vignettes deepSNV
	mv *.tar.gz builds
	
install: deepSNV.roxygen
	R CMD install deepSNV
	
clean: deepSNV*
	find deepSNV -name '.*' | xargs rm -rf
	find deepSNV -name '*.o' | xargs rm -rf
	find deepSNV -name '*.so' | xargs rm -rf
	rm deepSNV/config.*
	rm -rf autom4te.cache
	rm -rf deepSNV.Rcheck
	
bioc: roxygen clean
	#find  bioC/devel/deepSNV -type f | grep -v '/\.' | xargs rm 
	rsync -avp --exclude='.*' deepSNV/ bioc/devel/deepSNV/
	
