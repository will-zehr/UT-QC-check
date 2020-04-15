# UT-QC-check
Replicating QC steps to check if population stratification still happens

in directory */net/hunt/disk2/woverton/milewicz/burden_tests/new_QC*

- Step 1: include variants labeled PASS in vcf

```console
bcftools view -i 'FILTER="PASS"' /net/hunt/disk2/woverton/milewicz/burden_tests_withknowns/total_samples_combined.VQSR.final.vcf.gz > total_samples_combined.VQSR.final.passonly.vcf
bgzip total_samples_combined.VQSR.final.passonly.vcf
tabix -p vcf total_samples_combined.VQSR.final.passonly.vcf
```

- Step 2: Normalize

```console
bcftools norm -m- total_samples_combined.VQSR.final.passonly.vcf.gz -o total_samples_combined.passonly.VQSR.final.normalized.vcf
plink2 --vcf total_samples_combined.passonly.VQSR.final.normalized.vcf --out milewicz.VQSR --make-bed
```

- Step 3: filter out samples missing more than 10%

```console
vcftools --vcf total_samples_combined.passonly.VQSR.final.normalized.vcf --out missingness --missing-indv
```

This actually shows that no samples have more than 10% missingness.

sidestep: get rid of X and Y

```
awk 'BEGIN{FS=OFS="\t"} {gsub(/\X/, "23", $1)} 1' milewicz.VQSR.bim > test.bim
awk 'BEGIN{FS=OFS="\t"} {gsub(/\Y/, "24", $1)} 1' test.bim > test2.bim
mv test2.bim milewicz.VQSR.bim
```

and add sex/pheno information to .fam file

```python
import pandas as pd
import numpy as np
fam=pd.read_csv('milewicz.VQSR.fam', sep='\t', header=None)
fam=fam.rename(columns={1:'key'})
sex=pd.read_csv('../sample_ids/sample_sex.tsv', sep='\t', header=None)
sex=sex.rename(columns={0:'key'})
out=pd.merge(fam, sex, on='key')
out=out[[0,"key",2,3,1,5]]
out.isna().sum() #no NAN values
controls=pd.read_csv('../sample_ids/controls', header=None)
out['phenotype']=~out['key'].isin(controls[0]).astype(int)+3
out=out[[0,"key",2,3,1,'phenotype']]
out.to_csv('milewicz.VQSR.fam', sep='\t', header=False, index=False)
```

- Step 3: filter out very rare variants and those with more than 10% missingness

```console
plink --bfile milewicz.VQSR --maf  0.0000000000001 --out milewicz_passedQC_geno0.05_poly --geno 0.1 --make-bed >& plink_VQSR_call_rate_sample_filtering.log &
```

- Step 4: LD Prune

```console
plink --indep-pairwise 500 50 0.2 --bfile milewicz_passedQC_geno0.05_poly --out SAS_milewicz_passedQC_geno0.05_poly_500_50_0.2 >& plink_SAS_pruning.log &
```

and subset bfile with LD prune info

```
plink --bfile milewicz_passedQC_geno0.05_poly --extract SAS_milewicz_passedQC_geno0.05_poly_500_50_0.2.prune.in --make-bed --out milewicz_passedQC_geno0.05_poly_500_50_0.2.pruned >& plink_milewicz_pruning_2.log &
```

- Step 5: PCA

```console
plink2 --bfile milewicz_passedQC_geno0.05_poly_500_50_0.2.pruned --pca --out milewicz_geno0.05_pruned_allchrs &
```

## Checking for stratification


