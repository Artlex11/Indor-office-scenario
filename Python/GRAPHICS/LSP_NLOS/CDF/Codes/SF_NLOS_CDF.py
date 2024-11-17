import matplotlib.pyplot as plt
import numpy as np

SF_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\SF_NLOS_CDF.txt", 'r')

SF_NLOS = SF_NLOS.read()
numbers = list(map(float, SF_NLOS.split()))

# unique_numbers - массив уникальных значений
unique_numbers = []

# Формирование массива уникальных чисел

# unique_numbers = set(round(num, 7) for num in numbers)

for i in numbers:
    if i not in unique_numbers:
        unique_numbers.append(i)
    # else:
    #     print(i)

print(max(numbers))
print(min(numbers))

print(len(numbers))
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

sorted_unique_numbers = sorted(unique_numbers)
print(sorted_unique_numbers)

plt.plot(sorted(numbers), np.linspace(0,1, len(sorted(numbers))))
plt.title("CDF for SF_NLOS")
plt.xlabel("SF")
plt.ylabel("p")
plt.grid()
plt.show()
