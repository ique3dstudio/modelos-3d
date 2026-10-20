# Modelos 3D

Modelos em STL para imprimir em impressora 3D. Cada modelo tem um script gerador em `modelos/<nome>/` e os STLs prontos em `stl/<nome>/`.

## Balão giratório (`modelos/balao`)

Brinquedo sensorial para crianças com autismo (3 anos ou mais), inspirado no balão de losangos coloridos das fotos em `referencias/`. O balão tem 12 cm de diâmetro (150 mm de altura), fica pendurado na ponta arredondada de um pino e **gira livremente** quando a criança o empurra. Os losangos em relevo são para pintar com tinta acrílica não tóxica.

![Pré-visualização](modelos/balao/preview.png)

| Peça | Arquivo | Qtd | Tamanho (mm) | Observação |
|---|---|---|---|---|
| Envelope, metade A | `stl/balao/envelope_metade_A.stl` | 1 | 118 x 150 x 59 | oco, parede de 1,6 mm; já deitado na mesa, sem suporte |
| Envelope, metade B | `stl/balao/envelope_metade_B.stl` | 1 | 118 x 150 x 59 | idem; cole as duas metades |
| Cesta com colar | `stl/balao/cesta_com_colar.stl` | 1 | 35 x 35 x 62 | cesta trançada, 4 hastes grossas e colar que entra na boca |
| Pino | `stl/balao/pino.stl` | 1 | 172 x 10 x 10 | impresso deitado; precisa de mesa com pelo menos 180 mm |
| Base | `stl/balao/base.stl` | 1 | 120 x 120 x 30 | disco largo e pesado, para não tombar |

Não há peças pequenas soltas nem fios: todas as peças são grandes e coladas ou encaixadas.

### Montagem

1. Cole as duas metades do envelope (cola de plástico ou epóxi). Dentro, no topo, há um poste com uma covinha onde a ponta do pino apoia.
2. Cole o colar da cesta dentro da boca do envelope; o flange da cesta apoia na borda da boca.
3. Cole o pino no furo da base (ele entra 16 mm).
4. Abaixe o balão com a cesta sobre o pino, boca para baixo. A cesta passa em volta do colar da base e a ponta do pino encaixa na covinha. O balão gira apoiado só na ponta.

### Impressão

- PLA ou PETG, bico 0,4 mm, camada de 0,2 mm, 3 perímetros, sem suporte.
- Base com 100% de preenchimento (ou coloque peso dentro), para ficar firme.
- Pino com 100% de preenchimento e orientação deitada, como já vem no STL.
- Lixe as bordas e a ponta do pino, e confira a espessura de parede no fatiador.
- Os modelos foram checados no computador (peças fechadas, sem colisão entre elas, centro de massa abaixo do apoio), mas **não foram impressos nem testados**. Faça um teste de giro antes de dar a uma criança, e use sempre com a supervisão de um adulto: o brinquedo não é certificado.

### Regerar com outras medidas

```
pip install numpy scipy trimesh manifold3d networkx
python3 modelos/balao/gerar_balao.py --diametro 120 --parede 1.6
```
