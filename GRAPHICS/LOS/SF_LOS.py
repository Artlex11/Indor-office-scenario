import matplotlib.pyplot as plt

SF_LOS = open("D:\\C++\\Practice\\testing\\ConsoleApplication1\\LSP_for_LOS\\SF_LOS.txt", 'r')

SF_LOS = SF_LOS.read()
numbers = list(map(float, SF_LOS.split()))

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

# for i in unique_numbers:
#     grad.append(pow(10, i))

plt.scatter(unique_numbers, p)
plt.title("PDF for SF_LOS")
plt.xlabel("SF")
plt.ylabel("p")
plt.xlim(min(numbers), max(numbers))
plt.ylim(0, 1)
plt.grid()
plt.show()
