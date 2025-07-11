"""FHE Decryption and Result Interpretation for SecureGenomics Protocol."""

import tenseal as ts
from typing import Dict, Any

def decrypt_result(encrypted_result: bytes, private_crypto_context: bytes) -> Dict[str, Any]:
    private_crypto_context = ts.context_from(private_crypto_context)
    encrypted_result = ts.bfv_vector_from(private_crypto_context, encrypted_result)
    result = encrypted_result.decrypt()
    return result



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
]# same in encode.py

def interpret_result(result):
    # remove the extra one at the end to count the number of alleles
    result, num_genomes = result[:-1], result[-1]
    
    print('# Genomes:', num_genomes)
    
    num_alleles = num_genomes * 2
    allele_freqs = [r / num_alleles for r in result]
    
    allele_frequencies_map = {}
    for i, freq in enumerate(allele_freqs):
        allele_frequencies_map[f'{TARGET_VARIANTS[i][2]} ({TARGET_VARIANTS[i][0]})'] = {
            'odds_ratio': TARGET_VARIANTS[i][1],
            'allele_frequency': freq,
        }
    
    print('Allele frequencies:', allele_freqs)
    return {
        'num_genomes': num_genomes,
        'allele_frequencies_vector': allele_freqs,
        'allele_frequencies_map': allele_frequencies_map
    }
