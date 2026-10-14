#!/usr/bin/env python3
"""Gera os STLs do balão de ar quente (envelope em duas metades + cesta).

Uso:  python3 gerar_balao.py [--altura 120] [--parede 1.6] [--saida ../../stl/balao]

Dependências: numpy, scipy, trimesh, manifold3d  (pip install numpy scipy trimesh manifold3d)
Unidades: milímetros.
"""
import argparse
import os

import numpy as np
import trimesh
from scipy.interpolate import PchipInterpolator

# Perfil do envelope: (z/H, r/H), de baixo (boca) para cima (topo), medido na foto.
PERFIL = [
    (0.00, 0.100), (0.06, 0.120), (0.15, 0.180), (0.30, 0.300), (0.45, 0.385),
    (0.60, 0.400), (0.80, 0.3464), (0.906, 0.257), (0.9759, 0.1368),
    (0.9984, 0.0349), (1.00, 0.000),
]
N_LOSANGOS_VOLTA = 16   # losangos ao redor do balão
N_LOSANGOS_ALTURA = 6   # losangos da boca ao topo
RELEVO = 1.2            # mm de relevo dos losangos (para pintar)
N_Z, N_T = 180, 96      # resolução: pontos no perfil / segmentos por metade


def perfil_denso(H, n):
    z, r = np.array(PERFIL).T
    f = PchipInterpolator(z * H, r * H)
    zz = np.linspace(0, H, n)
    return zz, np.clip(f(zz), 0, None)


def deslocamento(theta, arco, comp_arco, relevo):
    """Relevo em losangos: 1 dentro do losango, 0 fora, com borda suave."""
    u = theta / (2 * np.pi) * N_LOSANGOS_VOLTA
    v = arco / comp_arco * N_LOSANGOS_ALTURA
    dd = np.abs((u % 1) - 0.5) + np.abs((v % 1) - 0.5)
    t = np.clip((0.42 - dd) / 0.06, 0, 1)
    return relevo * t * t * (3 - 2 * t)


def metade_envelope(H, parede, theta0):
    z, r = perfil_denso(H, N_Z)
    dz, dr = np.gradient(z), np.gradient(r)
    norm = np.hypot(dz, dr)
    nz, nr = dr / norm, -dz / norm            # normal externa (no plano r-z)
    # arco acumulado, para losangos de proporção parecida em qualquer altura
    arco = np.concatenate([[0], np.cumsum(np.hypot(np.diff(z), np.diff(r)))])
    zi = z - parede * nz
    ri = np.clip(r - parede * nr, 0, None)
    ri[-1] = 0
    zi[-1] = z[-1] - parede

    th = np.linspace(theta0, theta0 + np.pi, N_T + 1)
    T, A = np.meshgrid(th, arco)               # (N_Z, N_T+1)
    d = deslocamento(T, A, arco[-1], RELEVO)
    d[-8:] *= np.linspace(1, 0, 8)[:, None]    # relevo some perto do topo
    d[:3] *= np.linspace(0, 1, 3)[:, None]     # e na boca (borda lisa)
    Ro = r[:, None] + d * nr[:, None]
    Zo = z[:, None] + d * nz[:, None] + 0 * T
    Ri = np.broadcast_to(ri[:, None], T.shape)
    Zi = np.broadcast_to(zi[:, None], T.shape)

    def xyz(R, Z):
        return np.stack([R * np.cos(T), R * np.sin(T), Z], -1).reshape(-1, 3)

    nv = N_Z * (N_T + 1)
    V = np.vstack([xyz(Ro, Zo), xyz(Ri, Zi)])
    idx = lambda off, j, i: off + j * (N_T + 1) + i
    F = []
    for j in range(N_Z - 1):
        for i in range(N_T):
            a, b = idx(0, j, i), idx(0, j, i + 1)
            c, e = idx(0, j + 1, i + 1), idx(0, j + 1, i)
            F += [(a, b, c), (a, c, e)]                       # casca externa
            a, b = idx(nv, j, i), idx(nv, j, i + 1)
            c, e = idx(nv, j + 1, i + 1), idx(nv, j + 1, i)
            F += [(a, c, b), (a, e, c)]                       # casca interna
    for i in range(N_T):                                      # aro da boca
        a, b = idx(0, 0, i), idx(0, 0, i + 1)
        c, e = idx(nv, 0, i + 1), idx(nv, 0, i)
        F += [(a, c, b), (a, e, c)]
    for i in (0, N_T):                                        # faces de corte
        for j in range(N_Z - 1):
            a, b = idx(0, j, i), idx(0, j + 1, i)
            c, e = idx(nv, j + 1, i), idx(nv, j, i)
            F += [(a, b, c), (a, c, e)]
    m = trimesh.Trimesh(V, F, process=True)
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(m)
    return m


def furar_fio(m, H, parede, angulos, z_furo, diam=1.4):
    """Furos radiais na região da boca para amarrar os fios da cesta."""
    _, r = perfil_denso(H, N_Z)
    cil = []
    for a in angulos:
        c = trimesh.creation.cylinder(radius=diam / 2, height=H * 0.12, sections=24)
        c.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
        c.apply_translation([H * 0.12 / 2 + H * 0.06, 0, 0])
        c.apply_transform(trimesh.transformations.rotation_matrix(a, [0, 0, 1]))
        c.apply_translation([0, 0, z_furo])
        cil.append(c)
    return trimesh.boolean.difference([m] + cil, engine="manifold")


def deitar_na_mesa(m, lado):
    """Gira a metade para a face de corte ficar na mesa (sem suportes)."""
    s = 1 if lado == 0 else -1
    R = trimesh.transformations.rotation_matrix(s * np.pi / 2, [1, 0, 0])
    m.apply_transform(R)
    m.apply_translation(-m.bounds[0] * [0, 0, 1] - [m.centroid[0], m.centroid[1], 0])
    return m


def cesta(H):
    k = H / 120
    rb, rt, alt = 6.2 * k, 7.6 * k, 12 * k
    pts = [(0, 0), (rb, 0)]
    n = 8                                         # fileiras de trama
    for s in range(n):
        z0, z1 = alt * s / n, alt * (s + 1) / n
        r0 = rb + (rt - rb) * s / n
        r1 = rb + (rt - rb) * (s + 1) / n
        pts += [(r0 + 0.5 * k, z0 + (z1 - z0) * 0.25), (r1 + 0.5 * k, z0 + (z1 - z0) * 0.7),
                (r1, z1)]
    pts += [(rt - 1.4 * k, alt), (rt - 1.4 * k, 2.0 * k), (0, 2.0 * k)]   # copo interno
    m = trimesh.creation.revolve(np.array(pts), sections=96)
    furos = []
    for a in np.deg2rad([45, 135, 225, 315]):
        c = trimesh.creation.cylinder(radius=0.7, height=alt * 2, sections=20)
        c.apply_translation([5.9 * k * np.cos(a), 5.9 * k * np.sin(a), alt / 2])
        furos.append(c)
    return trimesh.boolean.difference([m] + furos, engine="manifold")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--altura", type=float, default=120.0, help="altura do envelope (mm)")
    ap.add_argument("--parede", type=float, default=1.6, help="espessura da parede (mm)")
    ap.add_argument("--saida", default=os.path.join(os.path.dirname(__file__), "..", "..", "stl", "balao"))
    a = ap.parse_args()
    os.makedirs(a.saida, exist_ok=True)
    H = a.altura
    pecas = {}
    for nome, lado in (("envelope_metade_A", 0), ("envelope_metade_B", 1)):
        m = metade_envelope(H, a.parede, np.pi * lado)
        ang = np.deg2rad([45, 135] if lado == 0 else [225, 315])
        m = furar_fio(m, H, a.parede, ang, z_furo=H * 0.05)
        pecas[nome] = deitar_na_mesa(m, lado)
    pecas["cesta"] = cesta(H)
    for nome, m in pecas.items():
        p = os.path.join(a.saida, f"{nome}.stl")
        m.export(p)
        x, y, z = m.extents
        print(f"{nome}: {x:.1f} x {y:.1f} x {z:.1f} mm  watertight={m.is_watertight}  "
              f"volume={m.volume/1000:.1f} cm3 -> {p}")


if __name__ == "__main__":
    main()
