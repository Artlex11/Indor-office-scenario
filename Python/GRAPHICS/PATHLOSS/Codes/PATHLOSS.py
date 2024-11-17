import matplotlib.pyplot as plt
import numpy as np

# график для LOS
PATHLOSS_LOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\LOS\\PATHLOSS_LOS.txt", "r")
d_LOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\LOS\\d_LOS.txt", "r")

# чтение файла со значениями PATHLOSS
PATHLOSS_LOS = PATHLOSS_LOS.read()
numbers_LOS = list(map(float, PATHLOSS_LOS.split()))

# чтение файла со значениями расстояний d
d_LOS = d_LOS.read()
numbers_d_LOS = list(map(float, d_LOS.split()))

print(len(numbers_d_LOS))
print(len(numbers_LOS))

plt.scatter(numbers_LOS, numbers_d_LOS)
plt.title("Pathloss for LOS")
plt.xlabel("d (m)")
plt.ylabel("Pathloss (dB)")
plt.grid()
plt.show()

# график для NLOS
PATHLOSS_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\NLOS\\PATHLOSS_NLOS.txt", "r")
d_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\NLOS\\d_NLOS.txt", "r")

PATHLOSS_NLOS = PATHLOSS_NLOS.read()
numbers_NLOS = list(map(float, PATHLOSS_NLOS.split()))

d_NLOS = d_NLOS.read()
numbers_d_NLOS = list(map(float, d_NLOS.split()))

print(len(numbers_d_NLOS))
print(len(numbers_NLOS))

# Построение совместного графика
plt.scatter(numbers_NLOS, numbers_d_NLOS)
plt.title("Pathloss for NLOS")
plt.xlabel("d (m)")
plt.ylabel("Pathloss (dB)")
plt.grid()
plt.show()
