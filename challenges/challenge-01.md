> Entangle two qubits. Each atom encodes a qubit in $|g\rangle$, $|r\rangle$, and the blockade is the entangling resource: it forbids $|rr\rangle$, so a global drive takes the product state $|gg\rangle$ to the maximally entangled Bell state $(|gr\rangle + |rg\rangle)/\sqrt{2}$ — the same interaction that underlies two-qubit gates on this hardware. Maximize the Bell-state fidelity at a strong interaction strength and a weak one, beating the reference square pulse at both.

**SYSTEM.** Two atoms at spacing $r$, both coupled to the same global drive. For $r < R_b$ the doubly excited state $|rr\rangle$ is shifted out of resonance by $C_6/r^6$ — the blockade.

**TASK.** Design $\Omega(t)$, $\delta(t)$ that prepare the Bell state $|\Phi\rangle = \left(|gr\rangle + |rg\rangle\right)/\sqrt{2}$ from $|gg\rangle$ at spacing $r_1$, then re-optimize at $r_2$.

**REGISTER.**

$$
r_1 = 5.0\ \mu\mathrm{m}, \qquad r_2 = 6.5\ \mu\mathrm{m}.
$$

**REFERENCE PULSE (THE BASELINE TO BEAT).** A resonant square pulse exploiting the collective $\sqrt{2}$-enhanced Rabi oscillation of the blockaded pair:

$$
\Omega(t) = \Omega_{\mathrm{ref}} = 6.283\ \mathrm{rad}/\mu\mathrm{s}\ (2\pi \times 1.0\ \mathrm{MHz}), \qquad \delta(t) = 0, \qquad t \in [0, T_{\mathrm{ref}}], \qquad T_{\mathrm{ref}} = \frac{\pi}{\sqrt{2}\Omega_{\mathrm{ref}}} \approx 354\ \mathrm{ns}
$$

(realized as 352 ns on the 4 ns clock). At $\Omega_{\mathrm{ref}}$ the blockade radius is $R_b \approx 7.2\ \mu\mathrm{m}$, so both spacings are blockaded — but not equally:

| Spacing | $V = C_6/\hbar r^6$ | $V/\Omega_{\mathrm{ref}}$ | Blockade quality |
|---|---|---|---|
| $r_1 = 5.0\ \mu\mathrm{m}$ | 55.4 rad/µs | 8.8 | strong — reference pulse is close to optimal |
| $r_2 = 6.5\ \mu\mathrm{m}$ | 11.5 rad/µs | 1.8 | weak — leakage to $\vert rr\rangle$; the gap you close |

**SCORE.** Bell-state fidelity $F = |\langle\Phi|\psi(T)\rangle|^2$ in simulation; hardware validation compares measured populations $P_{gg}, P_{gr}, P_{rg}, P_{rr}$ (recommended: 500 shots) against simulated values.

**SUCCESS.** $F$ above the reference-pulse value at **both** spacings, within the device envelope.

---

### Submission (Discord team channel by 16:30)

- [ ] Waveforms + register (Pulser sequence, JSON, or code that generates them)
- [ ] Simulated $F$ at both spacings, with the reference-pulse value you beat (same shots)
- [ ] Pasqal Cloud job IDs for your hardware validation runs
- [ ] Three sentences: what you changed relative to the baseline, and why it worked

Device limits: see the [device envelope](../README.md#device-envelope) — and **verify `Device.specs` at runtime**.
