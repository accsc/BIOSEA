# BIOSEA Repository

I've split the PubChem HTSFP in two files due the 100Mb limit per file in Github

To rebuild the original file:

- Uncompress both parts: 

bzip2 -d HTS_PUBCHEM_ZSCORE_Part01.tsv.bz2 HTS_PUBCHEM_ZSCORE_Part02.tsv.bz2

- Merge all in one file

cat HTS_PUBCHEM_ZSCORE_Part02.tsv >> HTS_PUBCHEM_ZSCORE_Part01.tsv
mv HTS_PUBCHEM_ZSCORE_Part01.tsv HTS_PUBCHEM_ZSCORE.tsv; rm HTS_PUBCHEM_ZSCORE_Part02.tsv

