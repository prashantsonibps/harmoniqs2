> The same idea at real scale. A recent benchmark ran exactly this MIS recipe on neutral-atom processors from ~10 to 85+ atoms and published how solution quality slides as systems grow. That curve is your opponent: rebuild the paper's instances at a size of your choosing and out-design its pulse schedule. Anything above the published line wins.

**REFERENCE.** [arXiv:2511.22967](https://arxiv.org/abs/2511.22967) benchmarks MIS solution quality for the quantum adiabatic algorithm and QAOA on random unit-disk graphs, on neutral-atom QPUs from ~10 to 85+ atoms, and reports solution quality degrading with system size.

**TASK.** Reproduce the paper's instances at one or more system sizes and optimize the full pulse schedule — waveform shapes, sweep profile, timing, and (where the instance permits) register geometry — to raise solution quality above the published curve, scored on the paper's own metric at matched instance size and shot count.

**BASELINE.** The published curve itself. As a starting schedule, the starter kit ships the Challenge 02 ramp stretched to $T = 6\ 000$ ns with the same $\Omega_b$, $\delta_0$, $\delta_f$; it will **not** beat the paper — it is scaffolding, not a strategy.

**SCORING METRIC.** Report the paper's solution-quality metric; where a self-contained number is needed, use the approximation ratio

$$
\mathcal{R} = \langle |S_{\mathrm{meas}}| \rangle / \alpha(G)
$$

over valid independent sets $S_{\mathrm{meas}}$, at matched shots.

**SUCCESS.** Solution quality above the published value at the same system size.

**BONUS.** Formulate a portfolio-optimization QUBO (background: [arXiv:2201.02773](https://arxiv.org/abs/2201.02773)), map it to an MIS instance embeddable as a unit-disk register, and run it through the same pipeline.

---

### Submission (Discord team channel by 16:30)

- [ ] Waveforms + register (Pulser sequence, JSON, or code that generates them)
- [ ] Simulated score (paper's metric / $\mathcal{R}$) with the published value you beat, matched instance size and shots
- [ ] Pasqal Cloud job IDs for your hardware validation runs
- [ ] Three sentences: what you changed relative to the baseline, and why it worked

Device limits: see the [device envelope](../README.md#device-envelope) — and **verify `Device.specs` at runtime**.
