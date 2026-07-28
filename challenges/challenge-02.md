> Scale from one entangled pair to a computation. Qubit placement programs the interaction graph (edge iff $r_{ij} < R_b$), the problem's solution is encoded in the ground state of the resulting many-body Hamiltonian, and an adiabatic sweep steers the register toward it; a computational-basis measurement returns the answer as a bitstring — a maximum independent set of the graph. Embed the target graph, then beat the linear ramp's success probability.

**SYSTEM.** Four to five atoms; positions are free parameters. The register induces a *unit-disk graph* $G$: vertices are atoms, and $(i,j)$ is an edge iff $r_{ij} < R_b$.

**TASK.** Choose positions so that $G$ equals a target graph below, then design a global sweep $\Omega(t)$, $\delta(t)$ from $|g \cdots g\rangle$ such that a final measurement returns a maximum independent set of $G$ — atoms found in $|r\rangle$ form the candidate set.

**TARGET GRAPHS.** Solve either (or both — score is per graph):

| Graph | Edges | $\alpha(G)$ | Embedding hint |
|---|---|---|---|
| $G_A$: star $K_{1,3}$ | $\lbrace (0,1), (0,2), (0,3)\rbrace $ | 3 | center atom plus three atoms at 120° and radius $\rho$, with $\rho < R_b < \sqrt{3}\rho$ |
| $G_B$: cycle $C_5$ | $\lbrace (0,1), (1,2), (2,3), (3,4), (4,0)\rbrace $ | 2 | regular pentagon of side $s$, with $s < R_b < 1.618\ s$ (the diagonal) |

**BASELINE (THE RAMP TO BEAT).** A linear detuning sweep at constant amplitude, on the starter-kit registers ($\rho = s = 5.5\ \mu\mathrm{m}$), total $T = 4\ 000$ ns:

$$
\Omega(t):\quad 0 \ \xrightarrow{\ 250\ \mathrm{ns\ rise}\ }\  \Omega_b \ \xrightarrow{\ 3\ 500\ \mathrm{ns\ hold}\ }\  \Omega_b \ \xrightarrow{\ 250\ \mathrm{ns\ fall}\ }\  0, \qquad \Omega_b = 6.283\ \mathrm{rad}/\mu\mathrm{s}\ (2\pi \times 1.0\ \mathrm{MHz}),
$$

$$
\delta(t):\quad \delta_0 = -12.57\ \mathrm{rad}/\mu\mathrm{s} \ \xrightarrow{\ \mathrm{linear\ over\ the}\ 3\ 500\ \mathrm{ns\ hold}\ }\  \delta_f = +12.57\ \mathrm{rad}/\mu\mathrm{s} \quad (\mp 2\pi \times 2.0\ \mathrm{MHz}),
$$

with $\delta$ held at $\delta_0$ during the rise and $\delta_f$ during the fall. Sanity check on the pentagon at $s = 5.5\ \mu\mathrm{m}$: nearest-neighbor interaction $U_{\mathrm{nn}} = 31.3\ \mathrm{rad}/\mu\mathrm{s} > \delta_f$, diagonal interaction $U_{\mathrm{diag}} = 1.7\ \mathrm{rad}/\mu\mathrm{s} < \delta_f$ — the final detuning sits inside the MIS window.

**SCORE.** $P_{\mathrm{MIS}}$: the probability that the measured Rydberg configuration is an independent set of $G$ of maximum size $\alpha(G)$, from the same shot budget (recommended: 500 shots).

**SUCCESS.** $P_{\mathrm{MIS}}$ strictly above the baseline ramp on the same graph.

---

### Submission (Discord team channel by 16:30)

- [ ] Waveforms + register (Pulser sequence, JSON, or code that generates them)
- [ ] Simulated $P_{\mathrm{MIS}}$ with the baseline-ramp value you beat, same shots
- [ ] Pasqal Cloud job IDs for your hardware validation runs
- [ ] Three sentences: what you changed relative to the baseline, and why it worked

Device limits: see the [device envelope](../README.md#device-envelope) — and **verify `Device.specs` at runtime**.
