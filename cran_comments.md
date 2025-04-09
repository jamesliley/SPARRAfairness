# Comments for submission of package SPARRAfairness to CRAN

This is the second version of this package. 

## R CMD check results

Tested on the following platforms

* Windows Server 2022, R-devel, 64 bit

Three NOTEs on Windows Server 2022, R-devel, 64 bit

```
Possibly misspelled words in DESCRIPTION:
  SPARRA (2:48, 36:18)
  counterfactual (36:524)
```

These words are spelled correctly

```
* checking top-level files ... NOTE
Non-standard file/directory found at top level:
  'cran_comments.md'
```

This NOTE is due to the presence of this comment file.


```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```


This NOTE is discussed in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503) and can be ignored.

No ERRORs or WARNINGs for any tests.

## Downstream dependencies

None; first release

## Other

The methods in this package were developed by the authors. We will soon submit a paper detailing the methods and place a link in the DESCRIPTION.

Thank you for your consideration of our package.
