VERSION := $(shell cat DESCRIPTION | awk '/Version/{print $$2}')

all: roxygen build

roxygen: 
	#R CMD roxygen -u deepSNV
	find man  | xargs rm -rf
	R -e 'library(roxygen2); roxygenize(".")'

	patch --no-backup-if-mismatch NAMESPACE patches/NAMESPACE.patch
	autoconf -o configure configure.ac
	
check: build
	R CMD check "deepSNV_$(VERSION).tar.gz"
	
build: roxygen clean
	R CMD build --no-build-vignettes .
	
install: build
	R CMD install "deepSNV_$(VERSION).tar.gz"
	
clean: 
	find . -name '*.o' | xargs rm -rf
	find . -name '*.so' | xargs rm -rf
	find . -name 'config.*' | xargs rm -rf
	find . -name 'symbols.rds' | xargs rm -rf
	find . -name 'deepSNV.Rcheck' | xargs rm -rf
	find . -name 'autom4te.cache' | xargs rm -rf
	
bioc: roxygen clean
	#find  bioC/devel/deepSNV -type f | grep -v '/\.' | xargs rm 
	#rsync -avp --exclude='.*' deepSNV/ bioc/devel/deepSNV/
	
