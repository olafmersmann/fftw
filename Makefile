.PHONEY: clean test check build install pkg

install: clean
	(cd pkg; autoconf; ./cleanup)
	R CMD INSTALL --no-multiarch --configure-args="--with-fftw3-prefix=/opt/local" pkg

test: install
	Rscript pkg/inst/unittests/runner.r

check: clean
	R CMD check pkg && rm -fR pkg.Rcheck

clean:
	(cd pkg; ./cleanup)
	rm -fR pkg/src/*.o pkg/src/*.so pkg.Rcheck .RData .Rhistory

pkg: clean
	echo "Date: $(date +%Y-%m-%d)" >> pkg/DESCRIPTION
	git log --no-merges -M --date=iso --format=medium pkg/ > pkg/ChangeLog
	R CMD build pkg
	R CMD build --binary pkg
	git checkout pkg/DESCRIPTION
	rm -f pkg/ChangeLog
