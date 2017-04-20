#!/usr/bin/python
# -*- coding: utf-8 -*-

from pylab import *


'''
seq........- sekwencja (string)
window....- rozmiar okna (int)
time........- slownik translacji kodonow specyficzny dla gatunku (dict)
path........- sciezka do zapisu pliku wyjsciowego
- powinna zawierac takze nazwe pliku wyjsciowego,
a nie tylko miejsce jego ulokowania (string)
'''

def plotv(
    seq,
    window,
    time,
    path,
    ):
    """rysuje wykres dla dowolnej sekwencji wprowadzonej przez uzytkownika,
	oraz dla dowolnej wielkosci okna. Wymaga podania specyficznego
	dla gatunku slownika czasow translacji kodonow."""

    # przygotuj sekwencje:

    i = 0
    orf = []
    while i < len(seq) - 2:
        orf.append(seq[i:i + 3])
        i += 3
    while orf[-1] == 'UAA' or orf[-1] == 'UGA' or orf[-1] == 'UAG':
        orf.pop()

    # policz sredni czas translacji w aa/sec

    E = 0.0
    for codon in orf:
        if time.has_key(codon):
            E = E + time[codon]
    me = 1000 / (E / len(orf))

    # przygotuj window

    if window > len(orf):
        window = len(orf) - 2

    # przygotuj os x

    if window == 1:
        ox = range(1, len(orf) + 1)
    else:
        ox = range(window / 2, len(orf) + 1 - window / 2)

    # przygotuj os y

    oy = []
    i = 0
    if window == 1:
        for codon in orf:
            if time.has_key(codon):
                oy.append(1000 / float(time[codon]))
            else:
                oy.append(0)
    else:
        hw = window / 2
        i = hw
        while i <= len(orf) - hw:
            seq = orf[i - hw:i + hw]
            y = 0
            for s in seq:
                if time.has_key(s):
                    y = y + float(time[s])  # najpierw dodajemy caly czas
            oy.append(1000 * window / y)
			'''potem przerabiamy go na predkosc,
			dzielac przebyty dystans przez ten czas.
			calosc mnozymy razy 1000 bo chcemy miec prednosc w aa/sec.
			'''
            i += 1

    # ustal szerokosc figure

    mox = len(orf)
    moy = max(oy)
    if mox < 500:
        width = max(mox / 10, 10)
    else:
        width = min(mox / 20, 50)

    # ustal czestosc znacznikow na osi ox

    r = len(str(mox))
    if r == 1:
        hop = 1
    elif r == 2:
        hop = 2
    elif r == 3:
        hop = 10
    elif r == 4:
        hop = 50
    else:
        hop = 100

    # rysuj figure

    if mox < 1000:
        figure(figsize=(width, max(moy / 10, 4.5)), dpi=100)
        xlabel('codon  (sliding window size = ' + str(window) + ')')
        ylabel('V [aa/sec]')
        xticks(range(0, mox, hop))
        title('Your sequence', ha='center')
        b1 = axhline(y=me, color='r', linewidth=1)
        if oy.count(0) != 0:
            k = oy.index(0)
            b2 = axvline(x=ox[k], color='green', linewidth=1,
                         linestyle='--')
            legend(
                [b1, b2],
                ['mean translation speed V', 'translational readthrough'
                 ],
                loc='upper left',
                bbox_to_anchor=(0, 1.05),
                ncol=2,
                fancybox=True,
                shadow=True,
                )
        else:
            legend(
                [b1],
                ['mean translation speed V'],
                loc='upper left',
                bbox_to_anchor=(0, 1.05),
                fancybox=True,
                shadow=True,
                )
        plot(ox, oy, linewidth=2, linestyle='-')
        xlim(0, mox + 1)
        ylim(0, moy + 4)
        grid(True)
    else:
        nr = mox / 1000 + 1
        figure(figsize=(width, max(nr * moy / 10, 4.5)), dpi=100)
        i = 1
        xp = 0
        yp = 0
        while i <= nr:
            subplot(nr, 1, i)
            if i == nr:
                xlabel('codon  (sliding window size = ' + str(window)
                       + ')')
            ylabel('V [aa/sec]')
            xticks(range(xp, xp + 1000, hop))
            suptitle('Your sequence')
            b1 = axhline(y=me, color='r', linewidth=1)
            if oy.count(0) != 0:
                k = oy.index(0)
                b2 = axvline(x=ox[k], color='green', linewidth=1,
                             linestyle='--')
            if i == 1 and oy.count(0) != 0:
                legend(
                    [b1, b2],
                    ['mean translation speed V',
                     'translational readthrough'],
                    loc='upper left',
                    bbox_to_anchor=(0, 1.05),
                    ncol=2,
                    fancybox=True,
                    shadow=True,
                    )
            elif i == 1 and oy.count(0) == 0:
                legend(
                    [b1],
                    ['mean translation speed V'],
                    loc='upper left',
                    bbox_to_anchor=(0, 1.05),
                    fancybox=True,
                    shadow=True,
                    )
            plot(ox[xp:xp + 1000], oy[yp:yp + 1000], linewidth=2)
            xlim(xp, xp + 1000)
            ylim(0, moy + 4)
            grid(True)
            xp = xp + 1000
            yp = yp + 1000
            i = i + 1
    savefig(path, dpi=100)
    close()
