program test_SF_SPIN
  USE SF_SPIN
  USE SF_LINALG, only:diag,kron
  USE ASSERTING
  implicit none

  complex(8),dimension(4,4) :: gamma03_,gamma00_
  complex(8),dimension(4,4) :: gamma03,gamma00

  gamma00_ = diag([1d0,1d0,1d0,1d0])
  gamma03_ = diag([1d0,-1d0,1d0,-1d0])

  gamma00 = kron(pauli_tau_0,pauli_sigma_0)
  gamma03 = kron(pauli_tau_0,pauli_sigma_3)

  call assert(gamma00,gamma00_,"KRON tau_0.sigma_0")
  call assert(gamma03,gamma03_,"KRON tau_0.sigma_3")

end program test_SF_SPIN
