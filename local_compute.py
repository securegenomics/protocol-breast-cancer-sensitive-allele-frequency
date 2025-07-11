"""VCF Encoding Module for SecureGenomics Protocol - Alzheimer's Disease Allele Frequency Analysis."""

from typing import List, Dict, Any, TextIO
import pysam


# These are the most clinically significant and population-variable SNPs for AD research
TARGET_VARIANTS_BREAST_CANCER = [
    ("rs2981579", 1.27),   # FGFR2 proxy SNP; OR≈1.27:contentReference[oaicite:20]{index=20}
    ("rs3803662", 1.28),   # TOX3 SNP; OR≈1.28
    ("rs1045485", 1.10),   # CASP8 variant; OR≈1.10
    ("rs2981582", 1.26),   # FGFR2 intron 2; OR≈1.26:contentReference[oaicite:21]{index=21}
    ("rs1801516", 1.10),   # ATM missense; OR≈1.10
    ("rs17879961", 1.50),  # CHEK2 I157T; OR≈1.50
    ("rs11515",    1.10),  # CDKN2A promoter SNP; OR≈1.10
    ("rs25487",    1.10),  # XRCC1 Arg399Gln; OR≈1.10
    ("rs80357713", 5.00),  # BRCA1 truncating (185delAG); OR large (~5+)
    ("rs11571833", 1.53),  # BRCA2 K3326*; OR≈1.53:contentReference[oaicite:22]{index=22}
    ("rs180177102", 3.00), # PALB2 truncating; OR≈3.0
    ("rs555607708", 2.34), # CHEK2 1100delC; OR≈2.34:contentReference[oaicite:23]{index=23}
] # ref https://chatgpt.com/s/dr_686c84c4b4b48191b3c779a6393f15f0

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

import numpy as np

def local_compute(vcf_path: str) -> List[int]:
    """
    Calculate PRS for a given VCF file.
    """
    encoded_data = encode_vcf(vcf_path)
    
    odds_ratios = [_[1] for _ in TARGET_VARIANTS]
    
    prs = np.dot(encoded_data, np.log(odds_ratios))
    
    return prs
    
    
    
    
    
    
    
    
