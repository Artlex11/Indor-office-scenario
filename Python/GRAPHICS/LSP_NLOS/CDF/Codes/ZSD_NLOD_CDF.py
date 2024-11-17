import matplotlib.pyplot as plt
import numpy as np

ZSD_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\ZSD_NLOS_CDF.txt", 'r')

ZSD_NLOS = ZSD_NLOS.read()
numbers = list(map(float, ZSD_NLOS.split()))

# unique_numbers - массив уникальных значений
unique_numbers = []

# Формирование массива уникальных чисел

for i in numbers:
    if i not in unique_numbers:
        unique_numbers.append(i)

print(max(numbers))
print(min(numbers))

print(len(unique_numbers))

# n - массив количества уникальных значений
# p - массив вероятностей уникальных чисел
n = []
p = []
# k - количество уникальных чисел
k = 0

for i in unique_numbers:
    for j in numbers:
        if i == j:
            k+=1
    n.append(k)
    p.append(k/len(numbers))
    k = 0

print(n)
print(p)
print(max(p))

grad = []

for i in unique_numbers:
    grad.append(pow(10, i))

plt.subplot(1, 2, 1)
plt.plot(sorted(unique_numbers), np.linspace(0,1, len(sorted(unique_numbers))))
plt.title("CDF for ZSD_NLOS (rad)")
plt.xlabel("Angle(rad)")
plt.ylabel("p")
# plt.xlim(min(numbers), max(numbers))
# plt.ylim(0, 1)
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(sorted(grad), np.linspace(0,1, len(sorted(grad))))
plt.title("CDF for ZSD_NLOS (deg)")
plt.xlabel("Angle(deg)")
plt.ylabel("p")
# plt.xlim(min(grad), max(grad))
# plt.ylim(0, 1)
plt.grid()

plt.show()
