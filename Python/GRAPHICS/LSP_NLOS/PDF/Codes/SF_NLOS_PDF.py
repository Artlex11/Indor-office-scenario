import matplotlib.pyplot as plt

SF_NLOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\SF_NLOS_PDF.txt", 'r')

SF_NLOS = SF_NLOS.read()
numbers = list(map(float, SF_NLOS.split()))

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

plt.scatter(unique_numbers, p)
plt.title("PDF for SF_NLOS")
plt.xlabel("SF")
plt.ylabel("p")
plt.xlim(min(numbers), max(numbers))
plt.ylim(0, 0.06)
plt.grid()
plt.show()
