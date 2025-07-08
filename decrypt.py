"""FHE Decryption and Result Interpretation for SecureGenomics Protocol."""

import tenseal as ts
from typing import Dict, Any

def decrypt_result(encrypted_result: bytes, private_crypto_context: bytes) -> Dict[str, Any]:
    private_crypto_context = ts.context_from(private_crypto_context)
    encrypted_result = ts.bfv_vector_from(private_crypto_context, encrypted_result)
    result = encrypted_result.decrypt()
    return result



TARGET_VARIANTS = [
    ("rs429358", 12.0, 'APOE ε4'),    # APOE ε4 homozygotes (late-onset AD OR ≈12:contentReference[oaicite:6]{index=6})
    ("rs7412", 0.62, 'APOE ε2'),      # APOE ε2 allele (protective; OR≈0.62, 38% reduced risk:contentReference[oaicite:7]{index=7})
    ("rs2075650", 4.18, 'TOMM40'),   # TOMM40 (tags APOE ε4); risk allele OR≈4.18:contentReference[oaicite:8]{index=8}
    ("rs199768005", 0.10, 'APOE rare'), # Rare APOE (p.V236E) variant; OR≈0.10 (highly protective):contentReference[oaicite:9]{index=9}
    ("rs6857", 1.27, 'NECTIN2'),      # NECTIN2 (APOE locus) variant; risk allele OR≈1.27
    ("rs11136000", 1.19, 'CLU'),  # CLU variant; OR≈1.19
    ("rs3851179", 1.18, 'PICALM'),   # PICALM variant; OR≈1.18
    ("rs6733839", 1.22, 'BIN1'),   # BIN1 variant; OR≈1.22
    ("rs6656401", 1.18, 'CR1'),   # CR1 variant; OR≈1.18
    ("rs3764650", 1.23, 'ABCA7'),   # ABCA7 variant; OR≈1.23
] # same in encode.py

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
