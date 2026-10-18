"""Box dependence of the one-band Loewdin dipoles vs the two-band Kane dipoles (N1/MS-02 ADQW 64/10/29 at 0 kV/cm and
the triple well 58/14/34/14/58 at 55 kV/cm; self-consistent, 5e11 cm^-2, 300 K). Three Loewdin variants:
  loc   - Loewdin over the localized states only (sp.prepare, used by N1/NL-02/NL-05/MS-02 up to qwkit 0.6.0)
  buf   - Loewdin over all solved states of the box (localized + box continuum, nsolve states), then keep localized
and the Kane two-component dipoles of the localized states (core.kane_dipoles). Writes box_scan.json."""
import json, numpy as np
from qwkit import sp, core as q

CASES = {"ADQW": ([(64, "well"), (10, "barrier"), (29, "well")], 0.0),
         "TW": ([(58, "well"), (14, "barrier"), (34, "well"), (14, "barrier"), (58, "well")], 55.0)}
out = {}
for name, (lay, F) in CASES.items():
    for spacer in (50.0, 100.0, 200.0, 300.0):
        r = sp.self_consistent(lay, F, ns_cm2=5e11, T=300.0, spacer=spacer, dz=0.5, doping="spacer", nstates=10, maxit=400)
        z, V, m, Eg, act = r["z"], r["V"], r["m"], r["Eg"], r["active"]
        dz = z[1] - z[0]; zc = 0.5 * (act[0] + act[1])
        for nsolve in (10, 24):
            E, psi = q.solve_states(z, V, m, Eg, nstates=nsolve)
            psi = psi / np.sqrt(np.sum(psi ** 2, axis=0) * dz)
            keep = q.localized(z, psi, act); n = min(int(keep.sum()), 4)
            ZL = q.dipoles(z - zc, q.orthonormalize(z, psi[:, keep])[:, :n])
            ZB = q.dipoles(z - zc, q.orthonormalize(z, psi)[:, keep][:, :n])
            ZK, ov = q.kane_dipoles(z, psi[:, keep][:, :n], E[keep][:n], V, m, Eg, zc)
            Ek = E[keep][:n]
            hw = np.arange(0.5, 1.5, 0.002) * (Ek[1] - Ek[0])
            pk = lambda Z: max(abs(q.chi2_shg(Ek, Z, h, 0.010, 1e18)) for h in hw)
            cL, cB, cK = pk(ZL), pk(ZB), pk(ZK)
            key = f"{name} spacer {spacer:.0f} nsolve {nsolve}"
            out[key] = dict(nloc=int(keep.sum()), z_loc=[ZL[0,1], ZL[1,2], ZL[0,2]], z_buf=[ZB[0,1], ZB[1,2], ZB[0,2]],
                            z_kane=[ZK[0,1], ZK[1,2], ZK[0,2]], chi_buf_over_loc=cB/cL, chi_kane_over_loc=cK/cL,
                            chi_kane_over_buf=cK/cB)
            o = out[key]
            print(f"{key:28s} nloc {o['nloc']} z13 loc {ZL[0,2]:.3f} buf {ZB[0,2]:.3f} kane {ZK[0,2]:.3f} | z23 {ZL[1,2]:.2f}/{ZB[1,2]:.2f}/{ZK[1,2]:.2f}"
                  f" | chi buf/loc {cB/cL:.4f} kane/loc {cK/cL:.4f} kane/buf {cK/cB:.4f}", flush=True)
json.dump(out, open("box_scan.json", "w"), indent=1, default=float)
