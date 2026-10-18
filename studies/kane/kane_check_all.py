"""Two-band Kane metric check for the III-V structures of N1, NL-02, NL-05 and MS-02. For each structure at zero bias
(self-consistent, 5e11 cm^-2, 300 K, 100 A barriers): largest one-band overlap of the raw nonparabolic states, largest
two-band overlap, dipoles z12, z23, z13 from Loewdin-orthonormalized localized states (sp.prepare, as in the studies)
and from the two-component Kane states (core.kane_dipoles), and the peak ground-state chi(2) (core.chi2_shg, gamma =
10 meV, N = 1e18 cm^-3 nominal) with both dipole sets."""
import numpy as np
from qwkit import sp, core as q

STRUCT = {"ADQW 58/14/34": [(58, "well"), (14, "barrier"), (34, "well")],
          "ADQW 64/10/29": [(64, "well"), (10, "barrier"), (29, "well")],
          "ADQW 64.5/18/36.5": [(64.5, "well"), (18, "barrier"), (36.5, "well")],
          "ADQW 70/12/32": [(70, "well"), (12, "barrier"), (32, "well")],
          "TW 58/14/34/14/58": [(58, "well"), (14, "barrier"), (34, "well"), (14, "barrier"), (58, "well")],
          "ADQW 50/15/30": [(50, "well"), (15, "barrier"), (30, "well")],
          "ADQW 62/20/29": [(62, "well"), (20, "barrier"), (29, "well")]}
for name, lay in STRUCT.items():
    r = sp.self_consistent(lay, 0.0, ns_cm2=5e11, T=300.0, spacer=100.0, dz=0.5, doping="spacer", nstates=10)
    z, psi, E, V, m, Eg, act = r["z"], r["psi"], r["E"], r["V"], r["m"], r["Eg"], r["active"]
    dz = z[1] - z[0]; zc = 0.5 * (act[0] + act[1])
    psi = psi / np.sqrt(np.sum(psi ** 2, axis=0) * dz)
    keep = q.localized(z, psi, act)
    S1 = psi[:, keep].T @ psi[:, keep] * dz
    Ek, Pk = E[keep][:4], psi[:, keep][:, :4]
    ZK, ov = q.kane_dipoles(z, Pk, Ek, V, m, Eg, zc)
    E4, P4, ZL, rho = sp.prepare(r, nkeep=4)
    hw = np.arange(0.5, 1.5, 0.002) * (E4[1] - E4[0])
    cL = max(abs(q.chi2_shg(E4, ZL, h, 0.010, 1e18)) for h in hw)
    cK = max(abs(q.chi2_shg(Ek, ZK, h, 0.010, 1e18)) for h in hw)
    pL = ZL[0, 1] * ZL[1, 2] * ZL[0, 2]; pK = ZK[0, 1] * ZK[1, 2] * ZK[0, 2]
    print(f"{name:20s} localized {keep.sum()} | overlap one-band {np.max(np.abs(S1 - np.diag(np.diag(S1)))):.3f} two-band {ov:.1e} | "
          f"z12 {ZL[0,1]:.2f}/{ZK[0,1]:.2f} z23 {ZL[1,2]:.2f}/{ZK[1,2]:.2f} z13 {ZL[0,2]:.2f}/{ZK[0,2]:.2f} | product ratio {pK/pL:.4f} | "
          f"peak chi2 ratio {cK/cL:.4f}", flush=True)
