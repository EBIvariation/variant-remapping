
- `reverse_strand.py`, included (uses [pysam](https://pysam.readthedocs.io/en/latest/api.html), allows reverse strand 
allele correction (see Note below for further explanation), and filtering based on bowtie2's alignment score)
- `replace_refs.py`, included (uses old REF allele in the case of new REF=ALT, these then get swapped (see Note 2))

**Note: reverse strand correction:**  
When a variant is mapped to the reverse strand, the corresponding allele in the output VCF file is reversed, aka 
converted to the forward strand, as alleles in VCF files are always described on the forward strand. For example:
INPUT VCF:
Old genome: `G (REF) > A (ALT)`
This maps onto the new genome on the reverse strand:
OUTPUT VCF:
New genome: `C (REF) > T (ALT)`

**Note 2: what happens when the new REF=ALT?**
Example: Original variant on the old reference sequence: `G (REF) > T (ALT)` 
Once remapped to the new reference sequence, we get: `T (REF) > T (ALT)`  
What happens is `replace_refs.py` will search for the old REF allele, and replace the new REF with it. So we get:  
`G (REF) > T (ALT)`  
And then `bcftools norm` will check to see if any REF alleles are different to the new reference genome, and if so, it 
will swap them:
`T (REF) > G (ALT)`  
This makes sense because we know that T is now the REF allele, which means the original G allele was technically a 
variant of this new reference genome.  
