import numpy as np
from VA import ppp_generate, proj_select_from_points, VA_RunLoop

lam, M, d, seed = 1.0, 5000, 10, 12345
pts, R = ppp_generate(lam, M, d, seed)
Aa0    = proj_select_from_points(pts)
res    = VA_RunLoop(pts, Aa0, order="max", max_iter=100000)

print("R_ball:", R)
print("Aa0   :", Aa0.tolist())
print("Aa*   :", res["Aa_final"].tolist())
print("iters :", res["iters"], "hit_inside:", res["hit_inside"])

