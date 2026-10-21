#!/usr/bin/env python3
"""Balão sensorial de gomos deslizantes (print-in-place) + suporte com aro.

O balão tem 16 gomos. Cada par de gomos vizinhos é ligado por uma junta de pino e olhal
que permite um deslize vertical de +-8 mm. Ao puxar um gomo para cima os outros ficam
pendurados em escada, como as folhas da foto de referência. Nenhum gomo se solta.
O gomo 15 e o gomo 0 não têm junta entre si (a corrente é aberta).

Uso:  python3 gerar_balao.py [--diametro 120] [--saida ../../stl/balao]

Dependências: numpy, scipy, trimesh, manifold3d, networkx.  Unidades: mm.
"""
import argparse
import os

import numpy as np
import trimesh
from scipy.interpolate import PchipInterpolator

# ----------------------------------------------------------------- parâmetros
N_GOMOS = 16
PAREDE = 2.0            # espessura da casca
GAP = 0.25              # folga de cada lado entre gomos (0,5 mm no total)
R_TOPO = 5.0            # raio da abertura no topo (os gomos terminam aqui)
RELEVO_MAX = 1.4        # maior altura de textura acima da superfície base
N_Z, N_T = 200, 30      # resolução da malha: linhas no perfil / colunas por gomo
LOSANGOS_ALTURA = 7     # losangos da boca ao topo
V0 = 0.5

# perfil (z/H, r/H): boca embaixo, topo em cima; diâmetro máximo = 0,8 H
PERFIL = [
    (0.00, 0.100), (0.06, 0.120), (0.15, 0.180), (0.30, 0.300), (0.45, 0.385),
    (0.60, 0.400), (0.80, 0.3464), (0.906, 0.257), (0.9759, 0.1368),
    (0.9984, 0.0349), (1.00, 0.000),
]

# junta pino + olhal
PINO_R, CABECA_R, CABECA_H = 2.2, 3.6, 1.5
OLHAL_FURO, OLHAL_R, OLHAL_H = 2.6, 4.6, 8.0
LIVRE = 24.0            # comprimento livre do pino entre as cabeças (curso = LIVRE - OLHAL_H)
PINO_Y = -4.4           # posição tangencial do pino (dentro do gomo k)
R_PINO = 44.0


def sm(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


# --------------------------------------------------------------------- perfil
class Perfil:
    def __init__(self, H):
        self.H = H
        z, r = np.array(PERFIL).T
        esc = (0.4 * H - RELEVO_MAX) / (0.4 * H)      # diâmetro final inclui o relevo
        f = PchipInterpolator(z * H, r * H * esc)
        zz = np.linspace(0, H, 6000)
        rr = np.clip(f(zz), 0, None)
        ok = rr >= R_TOPO
        zz, rr = zz[ok], rr[ok]
        s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(zz), np.diff(rr)))])
        self.comp = s[-1]
        sg = np.linspace(0, s[-1], N_Z)
        self.s = sg
        self.z = np.interp(sg, s, zz)
        self.r = np.interp(sg, s, rr)
        tr, tz = np.gradient(self.r, sg), np.gradient(self.z, sg)
        self.nr, self.nz = tz, -tr                    # normal externa (unitária)
        self.zi = self.z - PAREDE * self.nz
        self.ri = np.clip(self.r - PAREDE * self.nr, 0.5, None)

    def r_interno(self, z):
        return float(np.interp(z, self.zi, self.ri))


# ------------------------------------------------------------------- texturas
def textura(theta, arco, comp):
    """Altura (mm) da textura em cada ponto. Losangos escalonados como no logo."""
    u = theta / (2 * np.pi) * N_GOMOS
    v = arco / comp * LOSANGOS_ALTURA + V0
    p, q = u + v, u - v
    a, b = np.round(p / 2), np.round(q / 2)
    dp, dq = p - 2 * a, q - 2 * b
    dd = np.maximum(np.abs(dp), np.abs(dq))
    tipo = np.mod(np.floor((a + b) / 2) + (a - b), 4)

    moldura = 0.6 * sm((0.95 - dd) / 0.07)
    t0 = np.full_like(dd, 1.4)                                        # liso, alto
    t1 = 0.4 + 1.0 * np.clip(np.cos(2 * np.pi * 2.5 * dp) * np.cos(2 * np.pi * 2.5 * dq) * 1.6, 0, 1)
    t2 = 0.5 + 0.9 * (0.5 + 0.5 * np.cos(2 * np.pi * 3.0 * dp))      # listras
    t3 = np.full_like(dd, -0.6)                                       # afundado
    interno = np.select([tipo == 0, tipo == 1, tipo == 2], [t0, t1, t2], t3)
    m = sm((0.50 - dd) / 0.04)
    return moldura * (1 - m) + interno * m


# ---------------------------------------------------------------------- gomos
def casca_gomo(P, k):
    th0, th1 = k * 2 * np.pi / N_GOMOS, (k + 1) * 2 * np.pi / N_GOMOS
    s_o = th0 + np.arcsin(GAP / P.r)
    e_o = th1 - np.arcsin(GAP / P.r)
    s_i = th0 + np.arcsin(GAP / P.ri)
    e_i = th1 - np.arcsin(GAP / P.ri)
    f = np.linspace(0, 1, N_T + 1)[None, :]
    To = s_o[:, None] + (e_o - s_o)[:, None] * f
    Ti = s_i[:, None] + (e_i - s_i)[:, None] * f
    A = np.broadcast_to(P.s[:, None], To.shape)

    d = textura(To, A, P.comp)
    rampa = np.ones(N_Z)
    n = 8
    rampa[:n] = np.linspace(0, 1, n)
    rampa[-n:] = np.linspace(1, 0, n)
    d = d * rampa[:, None]
    Ro = P.r[:, None] + d * P.nr[:, None]
    Zo = P.z[:, None] + d * P.nz[:, None]
    Ri = np.broadcast_to(P.ri[:, None], Ti.shape)
    Zi = np.broadcast_to(P.zi[:, None], Ti.shape)

    def xyz(R, Z, T):
        return np.stack([R * np.cos(T), R * np.sin(T), Z], -1).reshape(-1, 3)

    nv = N_Z * (N_T + 1)
    V = np.vstack([xyz(Ro, Zo, To), xyz(Ri, Zi, Ti)])
    idx = lambda off, j, i: off + j * (N_T + 1) + i
    F = []
    for j in range(N_Z - 1):
        for i in range(N_T):
            a, b = idx(0, j, i), idx(0, j, i + 1)
            c, e = idx(0, j + 1, i + 1), idx(0, j + 1, i)
            F += [(a, b, c), (a, c, e)]
            a, b = idx(nv, j, i), idx(nv, j, i + 1)
            c, e = idx(nv, j + 1, i + 1), idx(nv, j + 1, i)
            F += [(a, c, b), (a, e, c)]
    for i in range(N_T):
        for j, inv in ((0, False), (N_Z - 1, True)):          # aros da boca e do topo
            a, b = idx(0, j, i), idx(0, j, i + 1)
            c, e = idx(nv, j, i + 1), idx(nv, j, i)
            F += [(a, b, c), (a, c, e)] if inv else [(a, c, b), (a, e, c)]
    for i in (0, N_T):                                         # faces laterais
        for j in range(N_Z - 1):
            a, b = idx(0, j, i), idx(0, j + 1, i)
            c, e = idx(nv, j + 1, i), idx(nv, j, i)
            F += [(a, b, c), (a, c, e)]
    m = trimesh.Trimesh(V, F, process=True)
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(m)
    return m


def hull(*ms):
    return trimesh.convex.convex_hull(np.vstack([m.vertices for m in ms]))


def cil(r, h, x=0, y=0, z=0, sec=32):
    c = trimesh.creation.cylinder(radius=r, height=h, sections=sec)
    c.apply_translation([x, y, z])
    return c


def caixa(sx, sy, sz, x, y, z):
    b = trimesh.creation.box(extents=[sx, sy, sz])
    b.apply_translation([x, y, z])
    return b


def junta(P, zc):
    """Peças da junta entre o gomo k (pino) e o gomo k+1 (olhal), no quadro da emenda:
    x = radial, y = tangencial (+ para o gomo k+1), emenda em y = 0."""
    za, zb = zc - LIVRE / 2, zc + LIVRE / 2

    def haste(z_ref, sentido):
        """Escora a 45 graus da cabeça até a parede do gomo k."""
        z = z_ref
        for _ in range(30):
            xw = P.r_interno(z)
            z = z_ref + sentido * (xw - R_PINO)
        xw = P.r_interno(z)
        pad = caixa(1.6, 6.0, 5.0, xw - 0.2, PINO_Y, z)
        return hull(cabeca_ref[z_ref], pad), pad

    cabeca_ref = {}
    cab_b = cil(CABECA_R, CABECA_H, R_PINO, PINO_Y, za - CABECA_H / 2)
    cab_t = cil(CABECA_R, CABECA_H, R_PINO, PINO_Y, zb + CABECA_H / 2)
    cabeca_ref[za - CABECA_H] = cab_b
    cabeca_ref[zb + CABECA_H] = cab_t
    esc_b, _ = haste(za - CABECA_H, -1)
    esc_t, _ = haste(zb + CABECA_H, +1)
    pino = cil(PINO_R, LIVRE + CABECA_H, R_PINO, PINO_Y, zc)   # termina dentro das cabeças
    partes_k = [pino, esc_b, esc_t]                           # as escoras já incluem as cabeças

    xw = P.r_interno(zc)
    anel = cil(OLHAL_R, OLHAL_H, R_PINO, PINO_Y, zc)
    pad = caixa(1.6, 6.0, OLHAL_H, xw - 0.2, 4.0, zc)
    olhal = hull(anel, pad)
    furo = cil(OLHAL_FURO, 30, R_PINO, PINO_Y, zc)
    olhal = trimesh.boolean.difference([olhal, furo], engine="manifold")
    return partes_k, olhal


def girar(m, ang):
    m = m.copy()
    m.apply_transform(trimesh.transformations.rotation_matrix(ang, [0, 0, 1]))
    return m


def montar_gomos(H):
    P = Perfil(H)
    zc = 0.6 * H
    pk, olhal = junta(P, zc)
    gomos = []
    for k in range(N_GOMOS):
        partes = [casca_gomo(P, k)]
        if k < N_GOMOS - 1:                                    # pino: emenda com o gomo k+1
            partes += [girar(m, (k + 1) * 2 * np.pi / N_GOMOS) for m in pk]
        if k > 0:                                              # olhal: emenda com o gomo k-1
            partes.append(girar(olhal, k * 2 * np.pi / N_GOMOS))
        g = trimesh.boolean.union(partes, engine="manifold") if len(partes) > 1 else partes[0]
        gomos.append(g)
        print(f"gomo {k:2d}: watertight={g.is_watertight} faces={len(g.faces)}", flush=True)
    return gomos, P


# -------------------------------------------------------------------- suporte
def suporte(P):
    """Base larga + cesta trançada + 4 hastes + aro com assento cônico para o balão."""
    base_r, base_h = 65.0, 8.0
    rb = 4.0
    pts = [(0, 0), (base_r - rb, 0)]
    for a in np.linspace(-np.pi / 2, 0, 8):
        pts.append((base_r - rb + rb * np.cos(a), rb + rb * np.sin(a)))
    pts += [(base_r, base_h - 1), (base_r - 1, base_h), (0, base_h)]
    base = trimesh.creation.revolve(np.array(pts), sections=128)

    alt, rb_c, rt_c = 15.0, 15.0, 17.0
    pts = [(0, base_h - 0.5), (rb_c, base_h - 0.5)]
    n = 8
    for s in range(n):
        z0, z1 = alt * s / n, alt * (s + 1) / n
        r0, r1 = rb_c + (rt_c - rb_c) * s / n, rb_c + (rt_c - rb_c) * (s + 1) / n
        pts += [(r0 + 0.6, base_h + z0 + (z1 - z0) * 0.25), (r1 + 0.6, base_h + z0 + (z1 - z0) * 0.7),
                (r1, base_h + z1)]
    pts += [(0, base_h + alt)]
    cesta = trimesh.creation.revolve(np.array(pts), sections=96)

    z_cesta = base_h + alt
    z_aro = z_cesta + 70.0
    r_boca = P.r[0]
    r_aro = 26.0
    # aro: face de cima plana, de baixo em 45 graus para imprimir sem suporte
    pts_aro = [(0, z_aro + 5), (r_aro, z_aro + 5), (r_aro, z_aro + 3), (r_aro - 5, z_aro - 2),
               (0, z_aro - 2)]
    aro = trimesh.creation.revolve(np.array(pts_aro), sections=96)
    # assento cônico dentro da boca do balão (folga de 0,8 mm)
    ri_boca = P.ri[0]
    pts_c = [(0, z_aro + 5), (ri_boca - 0.8 - 0.0, z_aro + 5), (ri_boca - 0.8 + 3.0, z_aro + 15),
             (0, z_aro + 15)]
    cone = trimesh.creation.revolve(np.array(pts_c), sections=96)

    hastes = []
    for a in np.deg2rad([45, 135, 225, 315]):
        p0 = np.array([12.0 * np.cos(a), 12.0 * np.sin(a), z_cesta - 1.0])
        p1 = np.array([(r_aro - 5.5) * np.cos(a), (r_aro - 5.5) * np.sin(a), z_aro])
        v = p1 - p0
        L = np.linalg.norm(v)
        h = trimesh.creation.cylinder(radius=2.6, height=L + 3, sections=24)
        zax = v / L
        T = trimesh.geometry.align_vectors([0, 0, 1], zax)
        h.apply_transform(T)
        h.apply_translation((p0 + p1) / 2)
        hastes.append(h)
    pecas = [base, cesta, aro, cone] + hastes
    for m in pecas:
        if m.volume < 0:
            m.invert()
    s = trimesh.boolean.union(pecas, engine="manifold")
    return s, z_aro + 5


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diametro", type=float, default=120.0, help="diâmetro máximo do balão (mm)")
    ap.add_argument("--saida", default=os.path.join(os.path.dirname(__file__), "..", "..", "stl", "balao"))
    a = ap.parse_args()
    os.makedirs(a.saida, exist_ok=True)
    H = a.diametro / 0.8
    gomos, P = montar_gomos(H)
    conj = trimesh.util.concatenate(gomos)
    conj.export(os.path.join(a.saida, "balao_16_gomos.stl"))
    print("balao_16_gomos:", conj.extents.round(1), "mm")
    s, z_apoio = suporte(P)
    s.export(os.path.join(a.saida, "suporte.stl"))
    print("suporte:", s.extents.round(1), "mm  watertight", s.is_watertight, " z do apoio:", round(z_apoio, 1))


if __name__ == "__main__":
    main()
