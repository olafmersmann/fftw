OSTYPE		:= $(shell uname -s)
GITVERSION	:= $(word 3,$(shell git --version))

## Darwin / OS X specific flags:
ifeq ($(OSTYPE),Darwin)
INSTALL_FLAGS := --install-args='--no-multiarch --configure-args="--with-fftw3-prefix=/opt/local"'
endif

.PHONEY: clean test check build install pkg

install: clean
	(cd pkg; autoconf; ./cleanup)
	echo $(GITVERSION)
	R CMD INSTALL $(INSTALL_FLAGS) pkg

test: install
	Rscript pkg/inst/unittests/runner.r

check: clean
	R CMD check $(INSTALL_FLAGS) pkg && rm -fR pkg.Rcheck

clean:
	(cd pkg; ./cleanup)
	rm -fR pkg/src/*.o pkg/src/*.so pkg.Rcheck .RData .Rhistory

pkg: clean
	echo "Date: $(date +%Y-%m-%d)" >> pkg/DESCRIPTION
	git log --no-merges -M --date=iso --pretty=medium pkg/ > pkg/ChangeLog
	R CMD build pkg
	git checkout pkg/DESCRIPTION
	rm -f pkg/ChangeLog
