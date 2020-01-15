## Test environments
* local OS X install, x86_64-apple-darwin15.6.0, R 3.6.2
* local ubuntu 18.04.3 LTS, x86_64-pc-linux-gnu, R 3.6.2
* ubuntu 18.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr
(https://github.com/wch/checkresults/blob/master/httr/r-release).
All packages that I could install passed except:

* Ecoengine: this appears to be a failure related to config on
  that machine. I couldn't reproduce it locally, and it doesn't
  seem to be related to changes in httr (the same problem exists
  with httr 0.4).
