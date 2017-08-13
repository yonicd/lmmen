## Test environments
* local OS X install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* I have run R CMD check on the NUMBER downstream dependencies.
 ==> devtools::check(args = c('--as-cran'))

Updating lmmen documentation
Loading lmmen
Loading required package: lmmlasso
Loading required package: emulator
Loading required package: mvtnorm
Loading required package: miscTools
Loading required package: penalized
Loading required package: survival
----------------------------------------------------------------------
This is a test release of the package 'lmmlasso'. If you have any questions or problems, do not hesitate to contact the author.
----------------------------------------------------------------------
Welcome to penalized. For extended examples, see vignette("penalized").
Setting env vars --------------------------------------------------------------
CFLAGS  : -Wall -pedantic
CXXFLAGS: -Wall -pedantic
Building lmmen ----------------------------------------------------------------
'/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file --no-environ  \
  --no-save --no-restore --quiet CMD build '/Users/jonathans/projects/lmmen'  \
  --no-resave-data --no-manual 

* checking for file ‘/Users/jonathans/projects/lmmen/DESCRIPTION’ ... OK
* preparing ‘lmmen’:
* checking DESCRIPTION meta-information ... OK
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
Removed empty directory ‘lmmen/.Rd2pdf52030’
* building ‘lmmen_1.0.tar.gz’

Setting env vars --------------------------------------------------------------
_R_CHECK_CRAN_INCOMING_ : FALSE
_R_CHECK_FORCE_SUGGESTS_: FALSE
Checking lmmen ----------------------------------------------------------------
'/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file --no-environ  \
  --no-save --no-restore --quiet CMD check  \
  '/var/folders/4_/xhs9__yd49l4v4j4wdg9f0wr0000gp/T//Rtmpvf9Xby/lmmen_1.0.tar.gz'  \
  --as-cran --timings --as-cran --no-manual 

* using log directory ‘/Users/jonathans/projects/lmmen.Rcheck’
* using R version 3.3.2 (2016-10-31)
* using platform: x86_64-apple-darwin13.4.0 (64-bit)
* using session charset: UTF-8
* using options ‘--no-manual --as-cran’
* checking for file ‘lmmen/DESCRIPTION’ ... OK
* this is package ‘lmmen’ version ‘1.0’
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package ‘lmmen’ can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking examples ... OK
Examples with CPU or elapsed time > 5s
       user system elapsed
lmmen 4.357   0.36   5.431
** found \donttest examples: check also with --run-donttest
* DONE
Status: OK



R CMD check results
0 errors | 0 warnings | 0 notes

R CMD check succeeded
 . 
  
* FAILURE SUMMARY

* All revdep maintainers were notified of the release on RELEASE DATE.
