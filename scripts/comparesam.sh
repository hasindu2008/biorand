
PREFIX=test/minimap_merge_na12878/
./biorand comparesam $PREFIX/unipart.sam $PREFIX/multipart_merge.sam > out.tsv
diff $PREFIX/prev/out.tsv out.tsv
diff $PREFIX/prev/mismatches.tsv mismatches.tsv
diff $PREFIX/prev/only_in_a.tsv only_in_a.tsv
diff $PREFIX/prev/only_in_b.tsv only_in_b.tsv