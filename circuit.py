"""FHE Computation Circuit for SecureGenomics Protocol."""

from typing import Dict, Any, List
import tenseal as ts

def compute(encrypted_datasets: List[bytes]) -> bytes:
    # multiplication depth: 1
    vectors = [ts.bfv_vector_from(context=None, data=_) for _ in encrypted_datasets]
    
    num_alleles = len(vectors) * 2
    
    # sum vectors and divide by num_alleles
    encrypted_result = sum(vectors) / num_alleles
    
    return encrypted_result.serialize()