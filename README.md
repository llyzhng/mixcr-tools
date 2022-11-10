# mixcr-tools

Helper scripts to process clone tables from [MiXCR](https://docs.milaboratories.com/)

Only dependency is R package `optparse`. 

1. Combine two clones tables (e.g. clones_IGL.tsv and clones_IGK.tsv)

     The readFraction and/or uniqueUMIFraction is recalculated for the merged table

```
Rscript combine-clone-tables.R -x $table1 -y $table2 -o $combinedTable
```

2. Collapse clones by VDJ and CDR3 nucleotide sequence 

     Collapsed clone table has columns: VDJ, nSeqCDR3, readCount (or uniqueUMICount), cloneFraction

```
Rscript collapse-clone-tables-by-VDJ-CDR3.R -i $infile -o $outfile
```
