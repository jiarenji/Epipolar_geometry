from scipy.spatial.transform import Rotation as R

def combine(pitch_a, yaw_a, roll_a, pitch_rel, yaw_rel, roll_rel=0.0):
    r_aircraft = R.from_euler('ZYX', [yaw_a, pitch_a, roll_a], degrees=True)
    r_sensor   = R.from_euler('ZYX', [yaw_rel, pitch_rel, roll_rel], degrees=True)
    r_total = r_aircraft * r_sensor
    yaw, pitch, roll = r_total.as_euler('ZYX', degrees=True)
    return pitch, yaw, roll, r_total.as_matrix()

# ---------- first set ----------
p1, y1, r1, M1 = combine(
    -12.451446469046, 194.007591357786, 10.111773714869,
    -26.5, 264.57
)

# ---------- second set ----------
p2, y2, r2, M2 = combine(
    -12.739436042960, 193.783031934525, 9.746181341861,
    -26.22, 264.76
)

print("First set (Pitch,Yaw,Roll):", p1, y1, r1)
print("Second set (Pitch,Yaw,Roll):", p2, y2, r2)