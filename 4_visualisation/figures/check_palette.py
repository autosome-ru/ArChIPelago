"""Pairwise CIEDE2000 between the 7 series colours, for normal vision and simulated deuteranopia, protanopia, tritanopia
(Machado et al. 2009, severity 1.0). Usage: python3 check_palette.py  (reads PALETTE from archi_style.py)."""
import itertools, os, numpy as np
try:
    from archi_style import PALETTE_ORDER, PALETTE
except ImportError:
    PALETTE_ORDER = PALETTE = None

M = {"deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
     "protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
     "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]}

def hex2lin(h):
    c = np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)

def lin2lab(c):
    m = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])
    xyz = m @ np.clip(c, 0, 1) / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])

def de2000(l1, l2):
    L1, a1, b1 = l1; L2, a2, b2 = l2
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2); Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360; h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dL, dC = L2 - L1, C2p - C1p
    dh = h2p - h1p
    if C1p * C2p == 0: dh = 0
    elif dh > 180: dh -= 360
    elif dh < -180: dh += 360
    dH = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dh / 2))
    Lb, Cbp = (L1 + L2) / 2, (C1p + C2p) / 2
    if C1p * C2p == 0: hb = h1p + h2p
    elif abs(h1p - h2p) <= 180: hb = (h1p + h2p) / 2
    else: hb = (h1p + h2p + 360) / 2 if h1p + h2p < 360 else (h1p + h2p - 360) / 2
    T = 1 - 0.17 * np.cos(np.radians(hb - 30)) + 0.24 * np.cos(np.radians(2 * hb)) + 0.32 * np.cos(np.radians(3 * hb + 6)) - 0.20 * np.cos(np.radians(4 * hb - 63))
    dth = 30 * np.exp(-((hb - 275) / 25) ** 2)
    Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    Sl = 1 + 0.015 * (Lb - 50) ** 2 / np.sqrt(20 + (Lb - 50) ** 2); Sc = 1 + 0.045 * Cbp; Sh = 1 + 0.015 * Cbp * T
    Rt = -np.sin(np.radians(2 * dth)) * Rc
    return float(np.sqrt((dL / Sl) ** 2 + (dC / Sc) ** 2 + (dH / Sh) ** 2 + Rt * (dC / Sc) * (dH / Sh)))

def report(order=PALETTE_ORDER, pal=PALETTE, verbose=True):
    worst = {}
    for vision in ["normal"] + list(M):
        labs = {}
        for k in order:
            c = hex2lin(pal[k])
            if vision != "normal": c = np.array(M[vision]) @ c
            labs[k] = lin2lab(c)
        pairs = sorted(((de2000(labs[a], labs[b]), a, b) for a, b in itertools.combinations(order, 2)))
        worst[vision] = pairs[0]
        if verbose: print(f"{vision:7s} min dE2000 = {pairs[0][0]:5.1f} ({pairs[0][1]} vs {pairs[0][2]});  next: " +
                          "; ".join(f"{a}/{b} {d:.1f}" for d, a, b in pairs[1:3]))
    return worst

if __name__ == "__main__":
    report()
    import re
    r = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "archi_style.R")).read()
    m = re.search(r'MODEL_COL <- setNames\(c\(([^)]*)\)', r)
    rcols = re.findall(r'#[0-9A-Fa-f]{6}', m.group(1))
    print("R and Python palettes identical:", [c.upper() for c in rcols] == [PALETTE[k].upper() for k in PALETTE_ORDER])
