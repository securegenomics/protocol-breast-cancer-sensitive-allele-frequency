"""VCF Encoding Module for SecureGenomics Protocol - Alzheimer's Disease Allele Frequency Analysis."""

from typing import List, Dict, Any, TextIO
import pysam

# Target variants for Alzheimer's disease sensitive allele frequency analysis
# These are the most clinically significant and population-variable SNPs for AD research
TARGET_VARIANTS = [
    ("rs429358", 12.0),    # APOE ε4 homozygotes (late-onset AD OR ≈12:contentReference[oaicite:6]{index=6})
    ("rs7412", 0.62),      # APOE ε2 allele (protective; OR≈0.62, 38% reduced risk:contentReference[oaicite:7]{index=7})
    ("rs2075650", 4.18),   # TOMM40 (tags APOE ε4); risk allele OR≈4.18:contentReference[oaicite:8]{index=8}
    ("rs199768005", 0.10), # Rare APOE (p.V236E) variant; OR≈0.10 (highly protective):contentReference[oaicite:9]{index=9}
    ("rs6857", 1.27),      # NECTIN2 (APOE locus) variant; risk allele OR≈1.27
    ("rs11136000", 1.19),  # CLU variant; OR≈1.19
    ("rs3851179", 1.18),   # PICALM variant; OR≈1.18
    ("rs6733839", 1.22),   # BIN1 variant; OR≈1.22
    ("rs6656401", 1.18),   # CR1 variant; OR≈1.18
    ("rs3764650", 1.23),   # ABCA7 variant; OR≈1.23
], # ref https://chatgpt.com/s/dr_686c84c4b4b48191b3c779a6393f15f0

def make_record_map(vcf_path):
    """Create mapping from variant identifiers to genotype values."""
    record_map = {}
    vcf_reader = pysam.VariantFile(vcf_path)
    for record in vcf_reader.fetch():
        # Get genotype value (0=ref/ref, 1=ref/alt, 2=alt/alt)
        # Handle missing genotypes gracefully
        try:
            gt = record.samples[0]['GT']
            if None in gt:
                alt_count = 0  # Treat missing as reference
            else:
                alt_count = sum(gt)
        except (KeyError, IndexError):
            alt_count = 0
            
        # Map by rsID if available
        if record.id and record.id != '.':
            record_map[record.id] = alt_count
            
        # Map by chromosomal coordinates
        if record.alts:
            key = (str(record.chrom), record.pos, record.alts[0])
            record_map[key] = alt_count
            
    return record_map

def encode_on_variant_list(record_map: dict, filter_list: list):
    """Encode specific variants from the target list."""
    for pos_or_id in filter_list:
        if isinstance(pos_or_id, tuple):
            # Genomic coordinate specification
            key = (str(pos_or_id[0]), pos_or_id[1], pos_or_id[2])
        elif isinstance(pos_or_id, str):
            # rsID specification
            key = pos_or_id
        else:
            raise ValueError(f"Invalid variant specification: {pos_or_id}")
            
        # Return genotype count (0, 1, or 2) or 0 if variant not found
        yield record_map.get(key, 0)

def encode_vcf(vcf_path: str) -> List[int]:
    """
    Encode VCF genomic data to integer vectors suitable for FHE.
    
    Returns list of integers representing genotype counts for each target variant:
    - 0: homozygous reference (no risk alleles)
    - 1: heterozygous (one risk allele) 
    - 2: homozygous alternate (two risk alleles)
    
    For privacy-sensitive allele frequency analysis of Alzheimer's disease variants.
    """
    record_map = make_record_map(vcf_path)
    variants_list = [_[0] for _ in TARGET_VARIANTS]
    encoded_data = list(encode_on_variant_list(record_map, variants_list))
    return encoded_data