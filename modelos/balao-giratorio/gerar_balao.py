#!/usr/bin/env python3
"""Gera os STLs do balão giratório (brinquedo sensorial).

Peças: envelope em duas metades (com mancal interno), cesta com colar, pino e base.
O envelope fica pendurado na ponta arredondada do pino e gira livremente.

Uso:  python3 gerar_balao.py [--diametro 120] [--parede 1.6] [--saida ../../stl/balao]

Dependências: numpy, scipy, trimesh, manifold3d, networkx
Unidades: milímetros.
"""
import argparse
import os

import numpy as np
import trimesh
from scipy.interpolate import PchipInterpolator

# Perfil do envelope: (z/H, r/H), de baixo (boca) para cima (topo), medido na foto.
# O diâmetro máximo é 0,8 H.
PERFIL = [
    (0.00, 0.100), (0.06, 0.120), (0.15, 0.180), (0.30, 0.300), (0.45, 0.385),
    (0.60, 0.400), (0.80, 0.3464), (0.906, 0.257), (0.9759, 0.1368),
    (0.9984, 0.0349), (1.00, 0.000),
]
N_LOSANGOS_VOLTA = 16   # losangos ao redor do balão
N_LOSANGOS_ALTURA = 6   # losangos da boca ao topo
RELEVO = 1.2            # mm de relevo dos losangos (para pintar)
N_Z, N_T = 180, 96      # resolução: pontos no perfil / segmentos por metade

# Mecanismo giratório
POSTE_R = 9.73          # raio do mancal (poste que desce do topo do envelope)
POSTE_Z = 0.74          # altura (z/H) da face de baixo do mancal
PINO_R = 5.0            # raio do pino
FOLGA = 0.25            # folga de encaixe (mm)
ALTURA_BOCA = 57.0      # altura da boca do envelope acima do topo da base
ENCAIXE_PINO = 16.0     # profundidade do pino dentro da base
BASE_D, BASE_H = 120.0, 12.0


def perfil_denso(H, n):
    z, r = np.array(PERFIL).T
    # o diâmetro máximo (2*0,4*H) já inclui o relevo dos losangos
    f = PchipInterpolator(z * H, r * H * (0.4 * H - RELEVO) / (0.4 * H))
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
    nz, nr = -dr / norm, dz / norm            # normal externa (no plano r-z)
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


def mancal(H, parede, theta0):
    """Poste cheio que desce do topo; tem uma covinha cônica onde a ponta do pino apoia."""
    z_base = POSTE_Z * H
    z_topo = H - 1.37                           # termina dentro da parede do topo
    poste = trimesh.creation.cylinder(radius=POSTE_R, height=z_topo - z_base, sections=64)
    poste.apply_translation([0, 0, (z_topo + z_base) / 2])
    cone = trimesh.creation.cone(radius=6.5, height=6.5, sections=48)   # ápice para cima
    cone.apply_translation([0, 0, z_base - 0.01])
    poste = trimesh.boolean.difference([poste, cone], engine="manifold")
    return poste


def deitar_na_mesa(m, lado):
    """Gira a metade para a face de corte ficar na mesa (sem suportes)."""
    s = 1 if lado == 0 else -1
    R = trimesh.transformations.rotation_matrix(s * np.pi / 2, [1, 0, 0])
    m.apply_transform(R)
    m.apply_translation(-m.bounds[0] * [0, 0, 1] - [m.centroid[0], m.centroid[1], 0])
    return m


def boca_raio(H):
    return perfil_denso(H, 2)[1][0]


def boca_raio_interno(H, parede):
    return boca_raio(H) - parede


def cesta_com_colar(H, parede):
    """Cesta (anel trançado) + 4 hastes grossas + flange + colar que entra na boca."""
    k = H / 150
    rb, rt, alt = 15.0, 17.0, 15.0
    furo = 12.5
    pts = [(furo, 0), (rb, 0)]
    n = 8                                         # fileiras de trama
    for s in range(n):
        z0, z1 = alt * s / n, alt * (s + 1) / n
        r0 = rb + (rt - rb) * s / n
        r1 = rb + (rt - rb) * (s + 1) / n
        pts += [(r0 + 0.6, z0 + (z1 - z0) * 0.25), (r1 + 0.6, z0 + (z1 - z0) * 0.7), (r1, z1)]
    pts += [(furo, alt), (furo, 0)]
    cesta = trimesh.creation.revolve(np.array(pts), sections=96)

    folga_hastes = 35.0                           # comprimento livre das hastes
    z_flange = alt + folga_hastes
    r_col = boca_raio_interno(H, parede) - 0.6    # entra na boca com folga p/ cola
    partes = [cesta]
    for a in np.deg2rad([45, 135, 225, 315]):
        h = trimesh.creation.cylinder(radius=1.6, height=folga_hastes + 4, sections=24)
        h.apply_translation([14.5 * np.cos(a), 14.5 * np.sin(a), alt - 2 + (folga_hastes + 4) / 2])
        partes.append(h)
    flange = trimesh.creation.annulus(r_min=9.2, r_max=boca_raio(H) + 1.5, height=2.0, sections=96)
    flange.apply_translation([0, 0, z_flange + 1.0])
    colar = trimesh.creation.annulus(r_min=9.2, r_max=r_col, height=10.0, sections=96)
    colar.apply_translation([0, 0, z_flange + 2.0 + 5.0])
    partes += [flange, colar]
    m = trimesh.boolean.union(partes, engine="manifold")
    return m, z_flange + 2.0                      # z da face de apoio na boca


def pino(H, z_ponta_acima_base):
    """Pino com ponta arredondada, deitado para impressão."""
    comp = z_ponta_acima_base - 1.0          # o pino apoia no fundo do furo, 1 mm acima da base
    r_ponta = 2.5
    haste = comp - 5.0 - r_ponta
    pts = [(0, 0), (PINO_R, 0), (PINO_R, haste), (r_ponta, haste + 5.0)]
    for ang in np.linspace(0, np.pi / 2, 12):
        pts.append((r_ponta * np.cos(ang), haste + 5.0 + r_ponta * np.sin(ang)))
    m = trimesh.creation.revolve(np.array(pts), sections=64)
    m.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
    m.apply_translation([-m.bounds[0][0], 0, PINO_R])
    return m, comp


def base():
    """Disco largo e baixo com bordas arredondadas e colar para o pino."""
    r, h = BASE_D / 2, BASE_H
    rb = 5.0                                      # raio da borda arredondada
    pts = [(0, 0), (r - rb, 0)]
    for ang in np.linspace(-np.pi / 2, 0, 10):
        pts.append((r - rb + rb * np.cos(ang), rb + rb * np.sin(ang)))
    for ang in np.linspace(0, np.pi / 2, 10):
        pts.append((r - rb + rb * np.cos(ang), h - rb + rb * np.sin(ang)))
    pts += [(14.0, h), (11.0, h + 3.0), (11.0, h + ENCAIXE_PINO + 2.0), (0, h + ENCAIXE_PINO + 2.0)]
    m = trimesh.creation.revolve(np.array(pts), sections=128)
    furo = trimesh.creation.cylinder(radius=PINO_R + FOLGA, height=ENCAIXE_PINO + 1.0, sections=64)
    furo.apply_translation([0, 0, h + ENCAIXE_PINO + 2.0 - (ENCAIXE_PINO + 1.0) / 2 + 0.01])
    return trimesh.boolean.difference([m, furo], engine="manifold")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diametro", type=float, default=120.0, help="diâmetro máximo do envelope (mm)")
    ap.add_argument("--parede", type=float, default=1.6, help="espessura da parede (mm)")
    ap.add_argument("--saida", default=os.path.join(os.path.dirname(__file__), "..", "..", "stl", "balao-giratorio"))
    a = ap.parse_args()
    os.makedirs(a.saida, exist_ok=True)
    H = a.diametro / 0.8
    pecas = {}
    for nome, lado in (("envelope_metade_A", 0), ("envelope_metade_B", 1)):
        theta0 = np.pi * lado
        casca = metade_envelope(H, a.parede, theta0)
        m = trimesh.boolean.union([casca, mancal(H, a.parede, theta0)], engine="manifold")
        lado_y = 1 if lado == 0 else -1                      # recorta na metade y>=0 ou y<=0
        corte = trimesh.creation.box(extents=[4 * H, 2 * H, 4 * H])
        corte.apply_translation([0, lado_y * H, H / 2])
        m = trimesh.boolean.intersection([m, corte], engine="manifold")
        pecas[nome] = deitar_na_mesa(m, lado)

    cesta, z_apoio = cesta_com_colar(H, a.parede)
    pecas["cesta_com_colar"] = cesta
    # altura da ponta do pino acima do topo da base
    z_ponta = ALTURA_BOCA + POSTE_Z * H + 5.5   # 1 mm abaixo do ápice da covinha
    pecas["pino"], comp_pino = pino(H, z_ponta)
    pecas["base"] = base()

    for nome, m in pecas.items():
        p = os.path.join(a.saida, f"{nome}.stl")
        m.export(p)
        x, y, z = m.extents
        print(f"{nome}: {x:.1f} x {y:.1f} x {z:.1f} mm  watertight={m.is_watertight}  "
              f"volume={m.volume/1000:.1f} cm3")
    print(f"altura da boca acima da base: {ALTURA_BOCA} mm; comprimento do pino: {comp_pino:.0f} mm")


if __name__ == "__main__":
    main()
