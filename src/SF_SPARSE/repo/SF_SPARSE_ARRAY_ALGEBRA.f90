! Compatibility facade for the implemented same-orientation products.
! Mixed CSR/CSC products need an explicit result-orientation API and are not provided.
MODULE SF_SPARSE_ARRAY_ALGEBRA
  USE SF_SPARSE, only: matmul
  implicit none
  private
  public :: matmul
END MODULE SF_SPARSE_ARRAY_ALGEBRA
