# Image-Processing
Image Processing using MPI

Paralelizare:
Procesul 0 citeste matricea de pixeli din fisier si face broadcast catre toate
celelalte procese. Trimitem matricea integrala pentru a avea acces usor la vecini.
Fiecare proces calculeaza ce puncte ii revin lui(impartind numarul de puncte
la nr de procese). Daca impartirea nu se face exact, ultimul proces calculeaza
surplusul de cateva puncte.
Fiecare proces va calcula noua valoarea a pixelilor repartizati lui aplicand
filtrul si va trimite inapoi catre 0 o matrice doar cu pixelii corespunzatori(
nu trimite si punctele care nu i-au fost repartizate).
Procesul 0 primeste de la fiecare proces matricele lor si le asambleaza
intr-o matrice de dimensiunea egala cu matricea initiala.

Daca sunt mai multe filtre, se reiau pasii de mai sus, doar ca broadcastul
se face cu rezultatul pasului anterior.

0 trimite intr-adevar toata matricea catre celelalte procese, dar restul
proceselor ii vor trimite lui 0 o matrice de dimensiuni mult mai mici,
eficientizand transferul de date.

Scalare pentru rorschach.pgm si landscape.pnm

rorschach.pgm:
Teste pentru aplicarea unui singur filtru(de exemplu smooth):

1 proces: 0.72 0.72 0.71 0.73 0.73 -> medie 0.72 s
2 procese: 0.58 0.61 0.59 0.59 0.58 -> medie 0.59 s
3 procese: 0.54 0.54 0.54 0.56 0.54 -> medie 0.54 s
4 procese: 0.53 0.53 0.54 0.53 0.54 -> medie 0.53 s

Teste pentru aplicarea mai multor filtre(bssembssem):
1 proces: 4.61 4.62 4.60 4.64 4.61 -> medie 4.61 s
2 procese: 2.88 2.89 2.95 2.86 2.86 -> medie 2.88 s
3 procese: 2.26 2.27 2.28 2.25 2.25 -> medie 2.26 s
4 procese: 2.05 2.02 2.01 2.00 1.98 -> medie 2.01 s

landscape.pnm
Teste pentru aplicarea unui singur filtru(de exemplu smooth):
1 proces: 1.81 1.82 1.83 1.87 1.83 -> medie 1.83 s
2 procese: 1.27 1.26 1.28 1.28 1.25 -> medie 1.26 s
3 procese: 1.04 1.08 1.07 1.07 1.06 -> medie 1.06 s
4 procese: 0.96 0.99 1.00 0.98 0.97 -> medie 0.98 s

Teste pentru aplicarea mai multor filtre(bssembssem):
1 proces: 15.66 15.65 15.63 15.91 15.72 -> medie 15.71 s 
2 procese: 9.12 9.17 9.12 9.18 9.13 -> medie 9.14 s
3 procese: 7.08 7.06 7.03 6.97 7.04 -> medie 7.03 s
4 procese: 6.11 6.12 6.14 6.11 6.18 -> medie 6.13 s

Am aplicat doar filtrul smooth pentru ca tipul filtrului nu influenteaza timpul de executie,
sunt aceleasi operatii, doar ca pe alte numere.
Pentru ca scalabilitatea sa fie mai evidenta pe imaginea alb negru cu un singur filtru 
am fi putut sa rulam program de un numar mare de ori pentru a creste numarul de operatii
(de exemplu sa aplicam filtrul de 1000 de ori).
