# Modelos 3D

Modelos em STL para imprimir em impressora 3D. Cada modelo tem um script gerador em `modelos/<nome>/` e os STLs prontos em `stl/<nome>/`.

## Balão de ar quente (`modelos/balao`)

Inspirado no balão de losangos coloridos (fotos em `referencias/`). Medidas padrão: envelope com 120 mm de altura e 96 mm de diâmetro.

| Peça | Arquivo | Qtd | Observação |
|---|---|---|---|
| Envelope, metade A | `stl/balao/envelope_metade_A.stl` | 1 | oco, parede de 1,6 mm; já deitado na mesa, sem suporte |
| Envelope, metade B | `stl/balao/envelope_metade_B.stl` | 1 | idem; cole as duas metades |
| Cesta | `stl/balao/cesta.stl` | 1 | 16 x 16 x 12 mm, com trama em relevo |

- O relevo de losangos (1,2 mm) é para pintar à mão com tinta acrílica.
- Há 4 furos de 1,4 mm perto da boca do envelope e 4 furos verticais de 1,4 mm na cesta. Passe linha nos furos, como no balão da foto.
- Sugestão de impressão: PLA, bico 0,4 mm, altura de camada 0,2 mm, 3 perímetros, sem suporte. Cesta com camada de 0,12 mm.
- Os modelos não foram impressos nem testados; confira a espessura de parede no fatiador.

Para outro tamanho ou parede, use o escala do fatiador ou regere:

```
pip install numpy scipy trimesh manifold3d networkx
python3 modelos/balao/gerar_balao.py --altura 150 --parede 2.0
```
