from sympy import symbols
from numpy import sqrt
(nx, ny, nz) = symbols('nx ny nz')
(Vsss, Vsps, Vpps, Vppp) = symbols('Vsss Vsps Vpps Vppp')
(Vsds, Vpds, Vpdp, Vdds, Vddp, Vddd) = symbols('Vsds Vpds Vpdp Vdds Vddp Vddd')

#[s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2]
HOPPING_INTEGRALS = [[None for _ in xrange(9)] for __ in xrange(9)]

HOPPING_INTEGRALS[0][0] = Vsss
HOPPING_INTEGRALS[0][1] = nx * Vsps
HOPPING_INTEGRALS[0][2] = ny * Vsps
HOPPING_INTEGRALS[0][3] = nz * Vsps
HOPPING_INTEGRALS[1][0] = - nx * Vsps
HOPPING_INTEGRALS[2][0] = - ny * Vsps
HOPPING_INTEGRALS[3][0] = - nz * Vsps
HOPPING_INTEGRALS[0][4] = sqrt(3) * nx * ny * Vsds
HOPPING_INTEGRALS[0][5] = sqrt(3) * ny * nz * Vsds
HOPPING_INTEGRALS[0][6] = sqrt(3) * nx * nz * Vsds
HOPPING_INTEGRALS[4][0] = sqrt(3) * nx * ny * Vsds
HOPPING_INTEGRALS[5][0] = sqrt(3) * ny * nz * Vsds
HOPPING_INTEGRALS[6][0] = sqrt(3) * nx * nz * Vsds
HOPPING_INTEGRALS[0][7] = sqrt(3) / 2. * (nx * nx - ny * ny) * Vsds
HOPPING_INTEGRALS[7][0] = sqrt(3) / 2. * (nx * nx - ny * ny) * Vsds
HOPPING_INTEGRALS[0][8] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vsds
HOPPING_INTEGRALS[8][0] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vsds

HOPPING_INTEGRALS[1][1] = nx * nx * Vpps + (1 - nx * nx) * Vppp
HOPPING_INTEGRALS[1][2] = - nx * ny * (Vppp - Vpps)
HOPPING_INTEGRALS[2][1] = - nx * ny * (Vppp - Vpps)
HOPPING_INTEGRALS[1][3] = - nx * nz * (Vppp - Vpps)
HOPPING_INTEGRALS[3][1] = - nx * nz * (Vppp - Vpps)
HOPPING_INTEGRALS[1][4] = sqrt(3) * nx * nx * ny * Vpds + (1 - 2 * nx * nx) * ny * Vpdp
HOPPING_INTEGRALS[1][5] = nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[1][6] = sqrt(3) * nx * nx * nz * Vpds + (1 - 2 * nx * nx) * nz * Vpdp
HOPPING_INTEGRALS[4][1] = - sqrt(3) * nx * nx * ny * Vpds - (1 - 2 * nx * nx) * ny * Vpdp
HOPPING_INTEGRALS[5][1] = nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[6][1] = - sqrt(3) * nx * nx * nz * Vpds - (1 - 2 * nx * nx) * nz * Vpdp
HOPPING_INTEGRALS[1][7] = sqrt(3) * nx * (nx * nx - ny * ny) * Vpds + nx * (1 - nx * nx + ny * ny) * Vpdp
HOPPING_INTEGRALS[7][1] = - sqrt(3) * nx * (nx * nx - ny * ny) * Vpds - nx * (1 - nx * nx + ny * ny) * Vpdp
HOPPING_INTEGRALS[1][8] = - sqrt(3) * nx * nz * nz * Vpdp - 0.5 * nx * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
HOPPING_INTEGRALS[8][1] = sqrt(3) * nx * nz * nz * Vpdp + 0.5 * nx * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

HOPPING_INTEGRALS[2][2] = ny * ny * Vpps + (1 - ny * ny) * Vppp
HOPPING_INTEGRALS[2][3] = - ny * nz * (Vppp - Vpps)
HOPPING_INTEGRALS[3][2] = - ny * nz * (Vppp - Vpps)
HOPPING_INTEGRALS[2][4] = sqrt(3) * ny * ny * nx * Vpds + (1 - 2 * ny * ny) * nx * Vpdp
HOPPING_INTEGRALS[4][2] = - sqrt(3) * ny * ny * nx * Vpds - (1 - 2 * ny * ny) * nx * Vpdp
HOPPING_INTEGRALS[2][5] = sqrt(3) * ny * ny * nz * Vpds + (1 - 2 * ny * ny) * nz * Vpdp
HOPPING_INTEGRALS[5][2] = - sqrt(3) * ny * ny * nz * Vpds - (1 - 2 * ny * ny) * nz * Vpdp
HOPPING_INTEGRALS[2][6] = nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[6][2] = - nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[2][7] = sqrt(3) / 2. * ny * (nx * nx - ny * ny) * Vpds - ny * (1 - ny * ny + nx * nx) * Vpdp
HOPPING_INTEGRALS[7][2] = - sqrt(3) / 2. * ny * (nx * nx - ny * ny) * Vpds + ny * (1 - ny * ny + nx * nx) * Vpdp
HOPPING_INTEGRALS[2][8] = - sqrt(3) * ny * nz * nz * Vpdp - 0.5 * ny * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
HOPPING_INTEGRALS[8][2] = sqrt(3) * ny * nz * nz * Vpdp + 0.5 * ny * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

HOPPING_INTEGRALS[3][3] = nz * nz * Vpps + (1 - nz * nz) * Vppp
HOPPING_INTEGRALS[3][4] = nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[4][3] = - nx * ny * nz * (sqrt(3) * Vpds - 2 * Vpdp)
HOPPING_INTEGRALS[3][5] = sqrt(3) * nz * nz * ny * Vpds + (1 - 2 * nz * nz) * ny * Vpdp
HOPPING_INTEGRALS[5][3] = - sqrt(3) * nz * nz * ny * Vpds - (1 - 2 * nz * nz) * ny * Vpdp
HOPPING_INTEGRALS[3][6] = sqrt(3) * nz * nz * nx * Vpds + (1 - 2 * nz * nz) * nx * Vpdp
HOPPING_INTEGRALS[6][3] = - sqrt(3) * nz * nz * nx * Vpds - (1 - 2 * nz * nz) * nx * Vpdp
HOPPING_INTEGRALS[3][7] = sqrt(3) / 2. * nz * (nx * nx - ny * ny) * Vpds - nz * (nx * nx - ny * ny) * Vpdp
HOPPING_INTEGRALS[7][3] = - sqrt(3) / 2. * nz * (nx * nx - ny * ny) * Vpds + nz * (nx * nx - ny * ny) * Vpdp
HOPPING_INTEGRALS[3][8] = sqrt(3) * nz * (nx * nx + ny * ny) * Vpdp - 0.5 * nz * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
HOPPING_INTEGRALS[8][3] = - sqrt(3) * nz * (nx * nx + ny * ny) * Vpdp + 0.5 * nz * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

HOPPING_INTEGRALS[4][4] = nx**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vddd
HOPPING_INTEGRALS[5][5] = nz**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nz**2 + ny**2) * Vddp + nx**2 * Vddd
HOPPING_INTEGRALS[6][6] = nx**2 * nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + nz**2) * Vddp + ny**2 * Vddd
HOPPING_INTEGRALS[4][5] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
HOPPING_INTEGRALS[5][4] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
HOPPING_INTEGRALS[4][6] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
HOPPING_INTEGRALS[6][4] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
HOPPING_INTEGRALS[5][6] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
HOPPING_INTEGRALS[6][5] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
HOPPING_INTEGRALS[6][7] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
HOPPING_INTEGRALS[7][6] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
HOPPING_INTEGRALS[5][7] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
HOPPING_INTEGRALS[7][5] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
HOPPING_INTEGRALS[4][7] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
HOPPING_INTEGRALS[7][4] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
HOPPING_INTEGRALS[4][8] = 0.5 * sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
HOPPING_INTEGRALS[5][8] = - 0.5 * sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
HOPPING_INTEGRALS[6][8] = - 0.5 * sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
HOPPING_INTEGRALS[8][4] = 0.5 * sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
HOPPING_INTEGRALS[8][5] = - 0.5 * sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
HOPPING_INTEGRALS[8][6] = - 0.5 * sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))

HOPPING_INTEGRALS[7][7] = 0.25 * (nx**2 - ny**2)**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vdds
HOPPING_INTEGRALS[8][8] = 0.75 * (nx**2 + ny**2)**2 * Vddd + 3 * (nx**2 + ny**2) * nz**2 * Vddp + 0.25 * (nx**2 + ny**2 - 2*nz**2)**2 * Vdds
HOPPING_INTEGRALS[7][8] = 0.25 * (nx**2 - ny**2) * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
HOPPING_INTEGRALS[8][7] = 0.25 * (nx**2 - ny**2) * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
