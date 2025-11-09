
import string
import numpy as np
from matplotlib import pyplot as diagram


#dane wejściowe
tablicaX = np.array ([-1, 0, 1, 2])
tablicaY = np.array ([5, 6, 4, 7])

#eliminacja Gaussa dla macierzy kwadratowych
def Gauss (A, b):
	n = len(A)
	for i in range(n-1):
		#szukamy mocnych przekątnych
		wiersz = i
		for j in range(i+1, n):
			if abs(A[j][i]) > abs(A[wiersz][i]):
				wiersz = j
		#zamianka
		if wiersz != i:
			b[i], b[wiersz] = b[wiersz], b[i]
			A[i], A[wiersz] = A[wiersz], A[i]
		#tango down! Eliminujemy dolny trójkącik
		for j in range(i+1, n):
			mnoznik = A[j][i] / A[i][i]
			for k in range(i+1, n):
				A[j][k] = A[j][k] - mnoznik * A[i][k]
			b[j] = b[j] - mnoznik * b[i]
	#rozwiązywańsko w górę
	x = [0] * n
	for i in range(n-1, -1, -1):
		x[i] = b[i]
		for j in range(i+1, n):
			x[i] = x[i] - A[i][j] * x[j]
		x[i] = x[i] / A[i][i]
	return np.flip(x)

def AproksymacjaNStopnia (tablicaX, tablicaY, N):
	#metoda 2-go stopnia ma mieć 3 rzędy i 3 kolumny, itd.
	A = np.empty([N+1,N+1])
	b = np.empty([N+1])
	#wypełnianie macierzy od tyłu (od najmniejszych potęg x do najwyższych)
	for i in range (N+1):
		k = i
		b[-i-1] = np.sum((tablicaX**k)*tablicaY)
		for j in range (N+1):
			A[-i-1][-j-1] = np.sum(tablicaX**k)
			k = k+1
	return (Gauss(A,b))

def AlgorytmHornera (ilorazy, tablicaX, x):
	dl = len(tablicaX) - 1
	y = ilorazy[dl]
	for i in range(1,dl+1):
		y = ilorazy[dl-i] + (x - tablicaX[dl-i]) * y
	return y

def InterpolacjaNewtona (tablicaX, tablicaY):
	dl=len(tablicaX)
	A = np.zeros ([dl,dl])
	for i in range (dl):
		A[i][0] = tablicaY[i]
	
	for i in range (1, dl):
		k=dl-i
		for j in range (k):
			A[j][i] = (A[j+1][i-1] - A[j][i-1]) / (tablicaX[i+j] - tablicaX[j])
	ilorazy = A[0]
	return (ilorazy)

def Ekstrapolacja (wspolczynnikiWielomianu, punktx):
	dl = len(wspolczynnikiWielomianu)
	wynik = 0
	for i in range (dl):
		wynik = wynik + punktx**i * wspolczynnikiWielomianu[i]
	return wynik

print ('Dane wejściowe:')
for i in range(len(tablicaX)):
	print (string.ascii_uppercase[i],': (',tablicaX[i],',',tablicaY[i],')')

ilorazy=InterpolacjaNewtona (tablicaX, tablicaY)
print ('Ilorazy różnicowe: ',ilorazy)

x=np.arange (tablicaX[0], tablicaX[-1]+0.01, 0.01)
y=AlgorytmHornera(ilorazy, tablicaX, x)

#współczynniki wielomianu
N=3
wspolczynnikiWielomianu = AproksymacjaNStopnia(tablicaX, tablicaY, N)

#wyświetlanie współczynników
print ('Wielomian interpolujący:')
print (np.polynomial.Polynomial(np.round(wspolczynnikiWielomianu,decimals=3)))

N=2
print ('Wielomian aproksymujący ',N,'-go stopnia:')
print(np.polynomial.Polynomial(np.round(AproksymacjaNStopnia(tablicaX, tablicaY, N),decimals=3)))

N=1
print ('Wielomian aproksymujący ',N,'-go stopnia:')
print(np.polynomial.Polynomial(np.round(AproksymacjaNStopnia(tablicaX, tablicaY, N),decimals=3)))

N=2
aproksymacja = AproksymacjaNStopnia(tablicaX, tablicaY, N)

#EKSTRAPOLACJA
punktx = -3
print ('Ekstrapolacja dla punktu: ', punktx, ' wynosi: ',np.round(Ekstrapolacja(wspolczynnikiWielomianu, punktx),decimals=3))

#konfiguracja diagramu
diagram.grid(linestyle='--') #styl siatki w tle
diagram.title ('Wielomian interpolacji Newtona') #tytuł
diagram.gca().set_aspect('equal') #zachowanie proporcji na osiach XY
diagram.plot(tablicaX, tablicaY,'r.',label='punkty wejściowe') #punkty, kolor czerwony, użycie kropki
diagram.plot(x, y, 'g', label='wielomian interpolacji')
diagram.plot(x, Ekstrapolacja(aproksymacja,x), 'b', label='wielomian aproksymujący')
diagram.legend(loc='upper left') #legenda w lewym górnym rogu
diagram.show()
