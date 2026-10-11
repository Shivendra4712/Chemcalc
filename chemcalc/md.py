"""2D Lennard-Jones molecular dynamics (reduced units: eps = sigma = m = k_B = 1)."""
from dataclasses import dataclass
import numpy as np


@dataclass
class MDResult:
    n: int
    box: float
    time: np.ndarray
    kinetic: np.ndarray
    potential: np.ndarray
    temperature: np.ndarray
    snapshots: list          # [(step, positions)]
    rdf_r: np.ndarray
    rdf_g: np.ndarray
    thermostat_steps: int
    dt: float

    @property
    def total(self):
        return self.kinetic + self.potential

    @property
    def energy_drift_per_particle(self) -> float:
        """Spread of total energy during the NVE (no thermostat) phase."""
        seg = self.total[self.thermostat_steps:]
        return float(np.ptp(seg) / self.n)


def _forces(pos, box, rc):
    d = pos[:, None, :] - pos[None, :, :]
    d -= box * np.round(d / box)
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    inv2 = 1.0 / r2
    inv6 = inv2 ** 3
    within = r2 < rc * rc
    shift = 4 * (rc ** -12 - rc ** -6)
    f_over_r = np.where(within, 24 * (2 * inv6 ** 2 - inv6) * inv2, 0.0)
    forces = (f_over_r[:, :, None] * d).sum(1)
    pe = 0.5 * np.where(within, 4 * (inv6 ** 2 - inv6) - shift, 0.0).sum()
    return forces, pe


def run_md(n_side=6, density=0.7, temperature=1.0, steps=1500, dt=0.004,
           rc=2.5, thermostat=True, seed=0, n_snapshots=3) -> MDResult:
    """Velocity-Verlet MD with periodic boundaries.
    First half: velocity-rescaling thermostat (equilibration). Second half: NVE."""
    n = n_side ** 2
    box = float(np.sqrt(n / density))
    rc = min(rc, box / 2)
    spacing = box / n_side
    gx, gy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    pos = (np.column_stack([gx.ravel(), gy.ravel()]) + 0.5) * spacing
    rng = np.random.default_rng(seed)
    vel = rng.normal(0, np.sqrt(temperature), (n, 2))
    vel -= vel.mean(0)
    dof = 2 * n - 2
    temp_of = lambda v: (v ** 2).sum() / dof
    vel *= np.sqrt(temperature / temp_of(vel))
    f, pe = _forces(pos, box, rc)

    ke_a, pe_a, t_a = [], [], []
    snap_steps = set(int(s) for s in np.linspace(0, steps, n_snapshots))
    snaps = []
    therm = steps // 2 if thermostat else 0
    nb = 60
    edges = np.linspace(0, box / 2, nb + 1)
    hist, frames = np.zeros(nb), 0
    iu = np.triu_indices(n, 1)

    for step in range(steps + 1):
        ke = 0.5 * (vel ** 2).sum()
        ke_a.append(ke); pe_a.append(pe); t_a.append(temp_of(vel))
        if step in snap_steps:
            snaps.append((step, pos.copy()))
        if step >= steps // 2 and step % 10 == 0:
            d = pos[:, None, :] - pos[None, :, :]
            d -= box * np.round(d / box)
            dist = np.sqrt((d ** 2).sum(-1))[iu]
            hist += np.histogram(dist, bins=edges)[0]
            frames += 1
        if step == steps:
            break
        vel += 0.5 * dt * f
        pos += dt * vel
        pos %= box
        f, pe = _forces(pos, box, rc)
        vel += 0.5 * dt * f
        if thermostat and step < therm and step % 10 == 0:
            vel *= np.sqrt(temperature / temp_of(vel))

    r = 0.5 * (edges[1:] + edges[:-1])
    dr = edges[1] - edges[0]
    g = 2 * hist / max(frames, 1) / (n * (n / box ** 2) * 2 * np.pi * r * dr)
    return MDResult(n, box, np.arange(steps + 1) * dt, np.array(ke_a), np.array(pe_a),
                    np.array(t_a), snaps, r, g, therm, dt)
