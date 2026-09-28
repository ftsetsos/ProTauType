# ProTauType

ProTauType assigns six-SNP **MAPT haplotype labels** from a phased VCF. It reports the allele on each chromosome, its six-base haplotype string, and the corresponding names in three published tables: Allen, Labbé, and Babić Leko.

The script looks up variants by **rsID in the VCF `ID` column**, so it does not depend on fixed genomic coordinates. It does require the following six SNPs, in this output order:

| Position in string | VCF ID |
| --- | --- |
| 1 | `rs1467967` |
| 2 | `rs242557` |
| 3 | `rs3785883` |
| 4 | `rs2471738` |
| 5 | `rs8070723` |
| 6 | `rs7521` |

## Requirements

- Python 3; the script uses only the standard library for ordinary `.vcf` and `.vcf.gz` input.
- A **phased, diploid**, single-sample or multi-sample VCF containing all six SNPs, with their rsIDs in `ID` and a complete `GT` for every sample at each SNP. Genotypes must use `|` (for example, `0|1`).
- `bcftools` for `.bcf` input or the optional `--region` argument. Region queries require an indexed input file.

## Run

```bash
python MAPT_three_study_genotyper.py phased_chr17.vcf.gz -o MAPT_haplotypes.tsv
```

For an indexed GRCh38 file, you can restrict the scan to a region containing the six SNPs (adjust the region and chromosome name for your VCF):

```bash
python MAPT_three_study_genotyper.py phased_chr17.vcf.gz \
  --region chr17:45900000-46050000 \
  -o MAPT_haplotypes.tsv
```

Run `python MAPT_three_study_genotyper.py --help` for the command options.

### Starting with unphased data

Phase a sufficiently broad chr17 region containing the MAPT locus with a suitable tool such as [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html), then pass its phased VCF to ProTauType. A schematic Beagle command is:

```bash
java -Xmx32g -jar beagle.jar gt=chr17_region.vcf.gz out=chr17_phased
python MAPT_three_study_genotyper.py chr17_phased.vcf.gz -o MAPT_haplotypes.tsv
```

Supply the genetic map and reference panel appropriate to your genome build and study when phasing. Phase with surrounding variants; six isolated SNPs provide little information for resolving phase. If the phasing step imputes missing genotypes, review those calls and retain the original variant and sample quality information: this script reads the final `GT`, not whether a call was observed or imputed.

## Output

The output is tab-separated, with **two rows per VCF sample**. `_A` and `_B` correspond to the first and second alleles in each phased `GT` (`0|1` means the reference allele is on `_A` and the alternate allele is on `_B` at that site).

| Column | Meaning |
| --- | --- |
| `sample_id` | Original VCF sample ID with `_A` or `_B` appended. |
| Six `rs...` columns | Called base on that chromosome at each SNP, in the order above. |
| `hap_as_bases` | Concatenation of the six bases in that order. |
| `allen` | Exact string's label in Allen et al., Table 2. |
| `labbe` | Exact string's label in Labbé et al., Table 3. |
| `leko` | Exact string's label in Babić Leko et al., Table 2. |
| `H1_or_H2` | `H1` for `A` at `rs8070723`, `H2` for `G` at `rs8070723`. |

For example, a phased `GGGCAA` chromosome is labeled `B (H1b)` / `H1B` / `H1B` in the three study columns and `H1` in the final column. An `AGGCGG` chromosome is `A (H2a)` / `H2` / `H2A` and `H2`. Labels vary between studies; they are retained as printed in each table.

An **empty study column** means that exact six-base string has no entry in that study's lookup table. It is not a failed genotype and does not imply that the chromosome does not exist. The `H1_or_H2` field is derived separately from `rs8070723`, even if study labels are empty.

## Input checks and interpretation

- All six rsID records must be present exactly once. An absent or repeated target ID stops the run. Check the VCF's `ID` field, build, and `--region` if this happens.
- Every sample needs a valid phased diploid `GT` at every target SNP. An unphased (`/`) or missing (`.`) genotype stops the run rather than producing a partial haplotype. Resolve missing calls during upstream QC or phasing if appropriate.
- If heterozygous target SNPs have conflicting, populated `PS` (phase set) values, the script stops. Where phase-set tags are missing, it cannot verify that all `|` genotypes share a consistent phase block.
- The script does not phase data, measure phasing confidence, directly assay the MAPT inversion, or assign a structural haplotype from sequence evidence. `H1_or_H2` is a call based on the `rs8070723` proxy allele.
- The Babić Leko table uses the MAPT intron 9 deletion (`del-In9`) as a marker. This implementation uses `rs8070723` in its six-base representation as a proxy; verify concordance with a directly typed deletion when that distinction matters.

The script validates the input before opening the output file. On an input error it exits with a message and does not write a partial results table.

## Sources

The lookup dictionaries transcribe allele-string labels from:

1. Allen et al. (2014), *Alzheimer's Research & Therapy*, Table 2. [DOI: 10.1186/alzrt268](https://doi.org/10.1186/alzrt268).
2. Labbé et al. (2016), *Parkinsonism & Related Disorders*, Table 3. [DOI: 10.1016/j.parkreldis.2016.06.010](https://doi.org/10.1016/j.parkreldis.2016.06.010).
3. Babić Leko et al. (2018), *Brain and Behavior*, Table 2. [DOI: 10.1002/brb3.1128](https://doi.org/10.1002/brb3.1128).

## License

MIT is the proposed license for this repository. Add a `LICENSE` file with the appropriate copyright holder before publishing.
