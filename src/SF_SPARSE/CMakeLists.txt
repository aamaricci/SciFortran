# Explicit real/complex implementations; MATRIX_SPARSE is the reference only.
set(SF_SPARSE_src
  SF_SPARSE_COMMON.f90
  SF_SPARSE_CSR_D.f90
  SF_SPARSE_CSR_Z.f90
  SF_SPARSE_CSC_D.f90
  SF_SPARSE_CSC_Z.f90
  SF_SPARSE.f90
)
add_library(SF_SPARSELIB OBJECT ${SF_SPARSE_src})
