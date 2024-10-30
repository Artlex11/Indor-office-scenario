#define _USE_MATH_DEFINES // для C++
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <set>

#include <math.h> 

class User {
public:
    int id;
    double x, y, z; // Координаты пользователя
    double theta;   // Зенитный угол
    double phi;     // Азимутальный угол
    double wavelength; // Длина волны
    double d;         // Расстояние между элементами антенны
    int numElementsX = 4; // 4 элемента в ширину
    int numElementsY = 2; // 2 элемента в длину 

    // Ориентация массива User 
    double bearing; // ΩUT,α
    double downtilt; // ΩUT,β 
    double slant; // ΩUT,γ 

    User(int id, double x, double y, double z, double lambda, double bearing, double downtilt, double slant) : id(id), x(x), y(y), z(z), wavelength(lambda), bearing(bearing), downtilt(downtilt), slant(slant) {
        d = lambda / 2.0; // Расстояние между элементами 
        calculateSphericalCoordinates();
    }

    // вызывается при создании user
    void calculateSphericalCoordinates() { // мб стоит переделать и сразу посчитать сферические для каждого и запомнить , будет логичнее 
        double r = std::sqrt(x * x + y * y + z * z);

        // theta = asin(z / r)
        theta = acos(z / r); // Зенитный угол 
        phi = atan2(y, x);   // Азимутальный угол
    }

    // выводится после создания user через цикл 2
    void printSphericalCoordinates() const {
        std::cout << "User  " << id << " - Theta: " << theta << ", Phi: " << phi << ", Bearing: " << bearing << ", Downtilt: " << downtilt << ", Slant: " << slant << std::endl;
    }


    // Метод для вычисления диаграмм направленности
    double antennaPattern(double thetaAngle, double phiAngle) const {  // нет кросс-поляризации , делать по Table 7.5-6 Part-2: Channel model parameters for RMa (up to 7GHz) and Indoor-Office , так же нужно сделать углы поворота антенны 
        double pattern = 0.0;
        double k = 2 * M_PI / wavelength; // Волновое число

        // Расчет ψ1 и ψ2
        double psi1 = k * d * sin(thetaAngle);
        double psi2 = k * d * sin(phiAngle);

        // Амплитудная диаграмма (можно адаптировать по необходимости)
        double Af = 1.0; // Пример амплитудной диаграммы

        // Расчет ДН для плоской решётки
        pattern = Af * (sin(numElementsX * psi1 / 2) / sin(psi1 / 2)) * (sin(numElementsY * psi2 / 2) / sin(psi2 / 2));

        return std::abs(pattern);
    }


    //поправка в расчёте углов
    // углы посчитаны неверно, distance = std::sqrt(pow(receiver.x - transmitter.x, 2) + pow(receiver.y - transmitter.y, 2) + pow(receiver.z - transmitter.z, 2));


    void calculateLOSAngles(const User& transmitter, const User& receiver, double& phiAOD, double& thetaZOD, double& phiAOA, double& thetaZOA) const {
        double dx = receiver.x - transmitter.x;
        double dy = receiver.y - transmitter.y;
        double dz = receiver.z - transmitter.z;
        double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Угол AOD (от передатчика к приемнику)
        thetaZOD = acos(dz / distance); // Зенитный угол 
        phiAOD = atan2(dy, dx);         // Азимутальный угол 

        // Угол AOA (от приемника к передатчику)
        thetaZOA = thetaZOD; // Зенитный угол 
        phiAOA = atan2(-dy, -dx); // Азимутальный угол
    }

};

class RicianChannel {
public:
    int numPaths;
    double roomWidth;
    double roomLength;
    double roomHeight;

    RicianChannel(int paths, double width, double length, double height)
        : numPaths(paths), roomWidth(width), roomLength(length), roomHeight(height) {} // Нужен к - фактор =7 , т.е. добавить шум ,т.е. AWGN ???

    std::complex<double> transmit(const User& transmitter, const User& receiver, std::complex<double> signal) {
        // Прямой компонент 
        double directDistance = calculateDistance(transmitter, receiver);
        std::complex<double> directComponent = std::polar(1.0 / directDistance, 0.0); // пока без вероятностей , а так  Table 7.4.2-1 LOS probability (Note: The LOS probability is derived with assuming antenna heights of 3m for indoor ??? у нас 1 метр ...)

        // Генерация случайных величин для отражённых путей 
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> amplitudeDist(0.0, 1.0); //пока так , но далее нужно через вероятности 
        std::uniform_real_distribution<> phaseDist(0.0, 2 * M_PI);

        std::complex<double> multipathComponent(0.0, 0.0);
        for (int i = 0; i < numPaths - 1; ++i) {
            double amplitude = amplitudeDist(gen);
            double phase = phaseDist(gen);
            multipathComponent += std::polar(amplitude, phase);
        }

        // Нормализация многолучевого компонента
        multipathComponent /= std::sqrt(numPaths - 1);

        std::complex<double> receivedSignal = (directComponent + multipathComponent) * signal;

        return receivedSignal;
    }

private:
    double calculateDistance(const User& a, const User& b) {
        return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
    }
};

//--------------------------------------------------------------------------------------------------------------------//
//-------------------------------ОСНОВНАЯ ПРОГРАММА-------------------------------------------------------------------//

int main() {
    // STEP 1 
    
    // Размеры комнаты
    double roomLength = 120.0; // Длина
    double roomWidth = 50.0;    // Ширина
    double roomHeight = 3.0;    // Высота 
    double wavelength = 0.3;     // Пример длины волны (в метрах)

    // Создание пользователей с случайными координатами 
    std::vector<User> users;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> xDist(0.0, roomLength);
    std::uniform_real_distribution<> yDist(0.0, roomWidth);
    double userHeight = 1.0; // Высота пользователей

    //--------ЦИКЛ 1--------- Создаёт user, генерирует bearing, downtilt и slant, создаёт все переменные и инициализирует некоторые из них, после возникают парные координаты и рассчёт сферических координат, добавлеие user в vector
    
    // Расчёт углов для перехода между системами координат // можно попробовать задать в классе user
    for (int i = 0; i < 12; ++i) {
        double bearing = (rand() % 360) * M_PI / 180.0; // Случайный угол направления (0-360 градусов)
        double downtilt = (rand() % 90) * M_PI / 180.0; // Случайный угол наклона (0-90 градусов)
        double slant = (rand() % 360) * M_PI / 180.0; // Случайный угол наклона (0-360 градусов)
        // объединить  юзерами
        users.emplace_back(i, xDist(gen), yDist(gen), userHeight, wavelength, bearing, downtilt, slant); // добавляет нового user в конец вектора
    }

    //--------ЦИКЛ 2---------
    // Вывод сферических координат пользователей
    for (const auto& user : users) {
        user.printSphericalCoordinates();
    }
    std::cout << std::endl;


    // Создание канала 
    RicianChannel channel(400, roomWidth, roomLength, roomHeight); // 400 количество лучей


    // Выбор шести случайных пар пользователей // исправленный (пары не повторяются)

    // Перемешка пользователей, чтобы индексы точно не повторялись, избавляемся от переидексации user в vector
    std::shuffle(users.begin(), users.end(), gen);

    std::uniform_int_distribution<> userDist(0, users.size() - 1); 
    std::vector<std::pair<int, int>> selectedPairs;

    for (int i = 0; i < 6; ++i) {
        int user1 = users[2 * i].id;
        int user2 = users[2 * i + 1].id;
        selectedPairs.emplace_back(user1, user2);
    }

    // Передача сигнала для каждой пары 
    std::complex<double> signal(1.0, 0.0); // Пример сигнала с амплитудой 1 и фазой 0 
    for (const auto& pair : selectedPairs) {
        const User& transmitter = users[pair.first];
        const User& receiver = users[pair.second];

        double phiAOD, thetaZOD, phiAOA, thetaZOA;
        transmitter.calculateLOSAngles(transmitter, receiver, phiAOD, thetaZOD, phiAOA, thetaZOA);

        std::complex<double> receivedSignal = channel.transmit(transmitter, receiver, signal);

        // Вычисление диаграммы направленности 
        double patternTransmitter = transmitter.antennaPattern(thetaZOD, phiAOD); //ДН передатчика
        double patternReceiver = receiver.antennaPattern(thetaZOA, phiAOA); // ДН приёмника

        std::cout << "Transmission from User " << transmitter.id << " to User " << receiver.id << " - Received Signal: " << receivedSignal
            << " | Transmitter Pattern: " << patternTransmitter << " | Receiver Pattern: " << patternReceiver << std::endl;

        std::cout << "LOS AOD: (" << phiAOD << ", " << thetaZOD << "), " << "LOS AOA: (" << phiAOA << ", " << thetaZOA << ")\n" << std::endl;

        // STEP 2
        // Вычисление вероятностей LOS
        double distance, P_LOS;
        distance = std::sqrt(pow(transmitter.x - receiver.x, 2) + pow(transmitter.y - receiver.y, 2));
        std::cout << "Distance between transmitter and receiver is " << distance << std::endl;
        if (distance <= 5)
            P_LOS = 1;
        else if (distance > 5 and distance <= 49)
            P_LOS = exp(-(distance - 5) / 70.8);
        else if (distance > 49)
            P_LOS = exp(-(distance - 49) / 211.7) * 0.54;
        std::cout << "LOS probability between users " << transmitter.id << " and " << receiver.id << " = " << P_LOS << std::endl << std::endl;

        // STEP 3
        // Вычисление Pathloss_LOS для каждой пары
        double path_loss_LOS;
        double nu = 1; // частота в ГГц
        path_loss_LOS = 32.4 + 17.3 * log10(distance) + 20 * log10(nu);

        std::cout << "PathLoss in LOS between users " << transmitter.id << " and " << receiver.id << " = " << path_loss_LOS << std::endl << std::endl;
    }

    return 0;
}

