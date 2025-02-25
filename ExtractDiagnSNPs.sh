# This script extracts SNP genotypes at 2 diagnostics SNP sites from a VCF file. The genotypes is represented by a number (0, 1 or 2) that represents the number of reference allele at that position.

# Define the VCF file
VCF_FILE="bd_intercept_ingroup_miss95_minDP6_Belgium.vcf.gz"

# Define the SNP locations (Notation: CHROM:POS, eg: LOCATIONS=("NC_064303.1:17298749" "NC_064303.1:23959938")
LOCATIONS=("NC_064303.1:17298749" "NC_064303.1:23959938")

# Extract sample IDs from the VCF file
SAMPLES=$(bcftools query -l "$VCF_FILE")

# Extract genotypes and encode as 2, 1, 0, or NA
bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r $(IFS=,; echo "${LOCATIONS[*]}") "$VCF_FILE" |
awk -v samples="$SAMPLES" '
 BEGIN {
    # Split sample IDs into an array
    split(samples, sample_ids, "\n");
    # Print the header
    printf "Sample ID\tSNP1\tSNP2\n";
}
{
    # Process each line (location)
    if (NR == 1) {
        # Store SNP1 genotypes
        for (i = 3; i <= NF; i++) {
            if ($i == "./.") {
                snp1[sample_ids[i-2]] = "NA";
            } else {
                split($i, gt, "/");
                if (gt[1] == "0" && gt[2] == "0") {
                    snp1[sample_ids[i-2]] = 2;
                } else if ((gt[1] == "0" && gt[2] == "1") || (gt[1] == "1" && gt[2] == "0")) {
                    snp1[sample_ids[i-2]] = 1;
                } else if (gt[1] == "1" && gt[2] == "1") {
                    snp1[sample_ids[i-2]] = 0;
                } else {
                    snp1[sample_ids[i-2]] = "NA";  # Handle unexpected genotypes
                }
            }
        }
    } else if (NR == 2) {
        # Store SNP2 genotypes and print the table
        for (i = 3; i <= NF; i++) {
            if ($i == "./.") {
                snp2[sample_ids[i-2]] = "NA";
            } else {
                split($i, gt, "/");
                if (gt[1] == "0" && gt[2] == "0") {
                    snp2[sample_ids[i-2]] = 2;
                } else if ((gt[1] == "0" && gt[2] == "1") || (gt[1] == "1" && gt[2] == "0")) {
                    snp2[sample_ids[i-2]] = 1;
                } else if (gt[1] == "1" && gt[2] == "1") {
                    snp2[sample_ids[i-2]] = 0;
                } else {
                    snp2[sample_ids[i-2]] = "NA";  # Handle unexpected genotypes
                }
            }
        }
        # Print the table rows
        for (sample in sample_ids) {
            printf "%s\t%s\t%s\n", sample_ids[sample], snp1[sample_ids[sample]], snp2[sample_ids[sample]];
        }
    }
}'