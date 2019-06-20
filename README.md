# biorand

My random bioinformatics tools.

## Building and usage

```sh
make
biorand [PROGRAM] [OPTIONS]
```


## Available programs and description

### filterfq

Filter a fastq file for reads with a sudden phred quality drop. Parameters are compile time flags at the moment. By default, all the reads containing a relative phred quality drop of 2.0 over a window of 200 bases (at least one dropped region) are output.

```
phred qual __________         __________
                     |_______|            
                     <------->
              qual drop region of ~200 bases         
```

```sh
Usage: biorand filterfq [OPTIONS] reads.fastq [reads2.fastq reads3.fastq ...]
OPTIONS:
        -t      Number of threads [12]
```

### filterpaf

Applies filters to a Pairwise mApping Format (PAF) file as specified through user parameters.

```sh
Usage: biorand filterpaf [OPTIONS] reads.paf reads.fastq

Options :
   -x STR              profile (martian/insert/del/trans)
   --qmin INT          minimum query gap
   --qmax INT          maximum query gap
   --tmin INT          minimum target gap
   --tmax INT          maximum target gap
   --gap-diff FLOAT    gap difference ratio (-1.0 to disable)
   --qual-thresh FLOAT relative phed quality drop in query (0.0 to disable)
   --w-size INT        window size outside the gap for relative phred score
   --bed STR           bed file for output
   --trans INT         threshold for transposons (-1 to disable)

Definitions :
                       tmin < target_gap < tmax
                                <-->
      target(ref)  : ----------------------------
                        |      |    |      |
      query(read)  :    --------------------
                                <-->
                        qmin < query_gap < qmax

gap difference ratio : output only if |query_gap-target_gap|<|target_gap*gap-diff|, ignore if gap-diff<0.0 [current gap-diff 0.1]
with qual drop (relative drop of qual-thresh, window outside gap w-size)
```

A few pre-set profiles are available to be specified through -x.

| profile     | description             | default parameters |
|-------------|-------------------------|--------------------|
| **martian** |  split mappings  (a similar gap in both target and query with the gapped sequence has a phred quality drop)       | --tmin 200 --tmax 5000 --qmin 200 qmax 5000  --gap-diff 0.1 --qual-thresh 2.0 --w-size 500 --trans -1.0 |
| **insert**  |  insertions             |  --tmin -200 --tmax 200 --qmin 200 qmax 5000  --gap-diff -1.0  --qual-thresh 0.0 --trans -1 |
| **del**     |  deletions              | tmin 200 --tmax 5000 --qmin -200 qmax 200  --gap-diff -1.0  --qual-thresh 0.0 --trans -1 |
| **trans**   |  insertions with a phred quality drop in the inserted sequence and the inserted sequence is mapped in the target to elsewhere (potential transposons)  | --tmin -200 --tmax 200 --qmin 200 qmax 5000  --gap-diff -1.0  --qual-thresh 0.0 --trans 500 |


Note : the readsIDs in *reads.fastq* and *reads.paf* should be in same order. If the *minimap2* aligner was executed on *reads.fastq*, the resultant *reads.paf* would have the reads in the same order.

### comparesam

Extensively compares two Sequence Alignment Map (SAM) files for differences in the alignment. Useful for comparing the alignments resultant from two different programs or the same program run with different parameters.  This tool was used for comparing NA12878 Nanopore reads mapped using *minimap2* with default parameters to that resulted from [minimap2 with partitioned indexes](https://doi.org/10.1038/s41598-019-40739-8).


```sh
Usage: biorand comparesam a.sam b.sam
```

Note : the samfiles *a* and *b* should have the reads in the same order and the multiple mappings for a given read ID should be adjacently located.


Let  *a* and *b* be two samfiles for the same set of reads,
but mapped with different mappers (or same mapper with different options).
comparesam will compare *a* and *b* and give mappings statistics such as the number of reads which are :
- unmapped in both
- correct - the read maps to the same location (locations overlaps if overlap based evaluation is set in the compile time)
- incorrect - the read DO NOT map to the same location (locations DO NOT overlap if overlap based evaluation is set)
- only mapped in  *a* (unique to *a*)
- only mapped in  *b* (unique to *b*)
- primary mapping in *a* is a supplementary mapping in *b*
- primary mapping in *b* is a secondary mappings in *b*

comparesam will also output (as tsv files and bed file) reads that: 1, mismatch between *a* and *b*; 2, unique to *a*; and 3, unique to *b*.

functionality is as follows :
- samfile *a* and *b* are sequentially read while loading all the mappings for a
particular read name at a time.
- For each loaded read, we compare the mappings between *a* and *b*. In *a* always the primary mapping is considered despite the value of
CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY compile time flags.
- In *b* supplementary and secondary mappings can also be considered if the compile time flags CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY are set (see comments in the [code](https://github.com/hasindu2008/biorand/blob/master/comparesam.c)).
- The comparison statistics are updated during the comparison and
the entries will be written to the tsv and bed file if required.
- Finally we print the comparison statistics.

### olp

Given an input fastq file, outputs pair-wise exact overlaps for all-vs-all reads

```sh
Usage: biorand olp in.fq out.paf read_len num_reads min_overlap max_overlap
```

- in.fq : input fastq
- out.paf : output [Pairwise mApping Format (PAF)](https://github.com/lh3/miniasm/blob/master/PAF.md) file
- read_len : the read length
- num_reads : the number of reads
- min_overlap : the minimum number of bases to be considered an overlap
- max_overlap : the maximum number of bases to be considered an overlap