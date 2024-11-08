import matplotlib.pyplot as plt

ASA_NLOS = open("D:\\C++\\Practice\\testing\\ConsoleApplication1\\LSP_for_NLOS\\ASA_NLOS.txt", 'r')

ASA_NLOS = ASA_NLOS.read()
numbers = list(map(float, ASA_NLOS.split()))

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

plt.scatter(unique_numbers, p)
plt.title("PDF for ASA_NLOS")
plt.xlabel("Angle(rad)")
plt.ylabel("p")
plt.xlim(min(numbers), max(numbers))
plt.ylim(0, 0.04)
plt.grid()
plt.show()

plt.scatter(grad, p)
plt.title("PDF for ASA_NLOS (deg)")
plt.xlabel("Angle(deg)")
plt.ylabel("p")
plt.xlim(min(grad), max(grad))
plt.ylim(0, 0.04)
plt.grid()
plt.show()