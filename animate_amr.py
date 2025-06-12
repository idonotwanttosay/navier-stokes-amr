import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import postprocess as pp

steps = pp.available_steps("rho")
if not steps:
    raise SystemExit("no output")

xs, ys = pp.read_grid("rho")
X, Y = np.meshgrid(xs, ys)

fig, ax = plt.subplots(figsize=(6,5))


def update(step):
    ax.clear()
    rho = pp.load_field("rho", step, xs, ys)
    c = ax.contourf(X, Y, rho, levels=40, cmap="viridis")
    # skip drawing AMR grid overlays
    ax.set_title(f"step {step}")
    ax.set_aspect("equal")
    return c.collections

ani = FuncAnimation(fig, update, frames=steps, blit=False)
ani.save("Result/animation.mp4", fps=10)
print("Saved Result/animation.mp4")
