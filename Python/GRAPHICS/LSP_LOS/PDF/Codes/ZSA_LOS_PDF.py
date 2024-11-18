import matplotlib.pyplot as plt

ZSA_LOS = open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\ZSA_LOS_PDF.txt", 'r')

ZSA_LOS = ZSA_LOS.read()
numbers = list(map(float, ZSA_LOS.split()))

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
plt.scatter(unique_numbers, p)
plt.title("PDF for ZSA_LOS")
plt.xlabel("Zenith(rad)")
plt.ylabel("p")
plt.xlim(min(numbers), max(numbers))
plt.ylim(0, 0.07)
plt.grid()

plt.subplot(1, 2, 2)
plt.scatter(grad, p)
plt.title("PDF for ZSA_LOS (deg)")
plt.xlabel("Zenith(deg)")
plt.ylabel("p")
plt.xlim(min(grad), max(grad))
plt.ylim(0, 0.07)
plt.grid()
plt.show()