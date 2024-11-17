import matplotlib.pyplot as plt
import numpy as np

# график для LOS
SF_LOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\SF_LOS_CDF.txt", 'r')

SF_LOS = SF_LOS.read()
numbers_LOS = list(map(float, SF_LOS.split()))

# unique_numbers - массив уникальных значений
unique_numbers_LOS = []

# Формирование массива уникальных чисел

# unique_numbers = set(round(num, 7) for num in numbers)

for i in numbers_LOS:
    if i not in unique_numbers_LOS:
        unique_numbers_LOS.append(i)
    # else:
    #     print(i)

print(max(numbers_LOS))
print(min(numbers_LOS))

print(len(numbers_LOS))
print(len(unique_numbers_LOS))

# n - массив количества уникальных значений
# p - массив вероятностей уникальных чисел
n_LOS = []
p_LOS = []
# k - количество уникальных чисел
k = 0

for i in unique_numbers_LOS:
    for j in numbers_LOS:
        if i == j:
            k += 1
    n_LOS.append(k)
    p_LOS.append(k/len(numbers_LOS))
    k = 0

print(n_LOS)
print(p_LOS)
print(max(p_LOS))

sorted_unique_numbers_LOS = sorted(unique_numbers_LOS)
print(sorted_unique_numbers_LOS)


# график для NLOS
SF_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\SF_NLOS_CDF.txt", 'r')

SF_NLOS = SF_NLOS.read()
numbers_NLOS = list(map(float, SF_NLOS.split()))

# unique_numbers - массив уникальных значений
unique_numbers_NLOS = []

# Формирование массива уникальных чисел
# unique_numbers = set(round(num, 7) for num in numbers)

for i in numbers_NLOS:
    if i not in unique_numbers_NLOS:
        unique_numbers_NLOS.append(i)
    # else:
    #     print(i)

print(max(numbers_NLOS))
print(min(numbers_NLOS))

print(len(numbers_NLOS))
print(len(unique_numbers_NLOS))

# n - массив количества уникальных значений
# p - массив вероятностей уникальных чисел
n_NLOS = []
p_NLOS = []
# k - количество уникальных чисел
k = 0

for i in unique_numbers_NLOS:
    for j in numbers_NLOS:
        if i == j:
            k += 1
    n_NLOS.append(k)
    p_NLOS.append(k/len(numbers_NLOS))
    k = 0

print(n_NLOS)
print(p_NLOS)
print(max(p_NLOS))

sorted_unique_numbers_NLOS = sorted(unique_numbers_NLOS)
print(sorted_unique_numbers_NLOS)


# Построение совместного графика
plt.plot(sorted(numbers_LOS), np.linspace(0, 1, len(sorted(numbers_LOS))), label="LOS")
plt.plot(sorted(numbers_NLOS), np.linspace(0, 1, len(sorted(numbers_NLOS))), label="NLOS")
plt.legend(loc='upper left')
plt.title("CDF for SF_LOS and SF_NLOS")
plt.xlabel("SF")
plt.ylabel("p")
plt.grid()
plt.show()
