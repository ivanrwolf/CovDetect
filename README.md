# CovDetect

###### B blocks detection tool

Retrieve regions with coverage above threshould between two different
experiments

**Requirements:**
* python v2.7
* pandas v0.20
* numpy v1.13

**usage:**
>python CovDetect.py [-h] [-bp BASEPAIRS] [-stdv STANDEV] input

**positional arguments:**

 * input
 >Input File

**optional arguments:**

* -h, --help
 
>Show help message and exit

* -bp BASEPAIRS, --basepairs BASEPAIRS
>Determine base pairs for sliding bases

* -stdv STANDEV, --standev STANDEV
>Number of standard deviations used to select the bases
