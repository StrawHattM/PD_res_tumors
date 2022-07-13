For this Combined Analysis, I have esentially merged the .cns (segment) files that were output by Genevia into combined segment files.

.cns files:

 - common: Segments that are common between files have been taken like one single segment. Segments present in only 1 or 2 samples are still considered.
 - union: Segments that are common between files have been taken like multiple segments. Segments present in only 1 or 2 samples are still considered.
 - noCTRL: Segments have been filtered (by matching chromosome, start and end ) using D2 and D3 controls.
 - Exclusive: Segments that are exclusively in either CIS or CARBO (by matching chromosome, start and end). This have been generated using UNION cns files.

gene reports:

 - m2: two probes (copy number variated segments) needed for the gene to be reported
 - m3: three probes (copy number variated segments) needed for the gene to be reported

This is relevant for filtering and "hit" selection. The default filtering to m3 accounts for "sequencing wise confidence".
If we want to further apply a filter based on log2/strict integer copy number, it might be better to use the m2 file to start from a less-filtered dataset.
