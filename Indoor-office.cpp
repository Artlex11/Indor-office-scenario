#define _USE_MATH_DEFINES // Для C++
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <set>
#include <algorithm> // Для std::sort

// Функция распределения Лапласа
double laplaceDistribution(double x, double mu, double b) {
    if (b == 0) {
        throw std::invalid_argument("Parameter b must be non-zero.");
    }
    return (0.5 / b) * exp(-std::abs(x - mu) / b);
}

// Функция для генерации случайных значений x по нормальному распределению
double generateNormalRandom(double mean, double stddev, std::mt19937& gen) {
    std::normal_distribution<double> dis(mean, stddev);
    return dis(gen); // Генерируем случайное значение по нормальному распределению
}

class LargeScaleParameters {
public:
    double shadowFading; // Затенение 
    double riceanK; // Фактор Райса 
    double delaySpread; // Задержка распространения 
    double azimuthSpreadDeparture; // Угловой спред выхода
    double azimuthSpreadArrival; // Угловой спред входа 
    double zenithSpreadDeparture; // Зенитный спред выхода 
    double zenithSpreadArrival; // Зенитный спред входа 
    LargeScaleParameters(double sf, double k, double ds, double asd, double asa, double zsd, double zsa)
        : shadowFading(sf), riceanK(k), delaySpread(ds), azimuthSpreadDeparture(asd),
        azimuthSpreadArrival(asa), zenithSpreadDeparture(zsd), zenithSpreadArrival(zsa) {}
};

class UserTerminal {
public:
    int id;
    double x, y, z; // Координаты пользователя 
    double wavelength; // Длина волны 
    double d; // Параметры антенны 
    int numElementsX = 4; // 4 элемента антенны по X
    int numElementsY = 2; // 2 элемента антенны по Y 
    double bearingAngle; // Угол поворота 
    double downtiltAngle; // Угол наклона 
    double slantAngle; // Угол наклона

    

    UserTerminal(int userId, double xCoord, double yCoord, double zCoord, double lambda, double bearing, double downtilt, double slant)
        : id(userId), x(xCoord), y(yCoord), z(zCoord), wavelength(lambda), bearingAngle(bearing), downtiltAngle(downtilt), slantAngle(slant) {
        d = lambda / 2.0; // Расстояние между элементами антенны 
    }

    // Метод для вычисления ДН
    double antennaPattern(double thetaAngle, double phiAngle) const {
        double pattern = 0.0;
        double k = 2 * M_PI / wavelength; // Волновое число

        // Расчет psi1 и psi2 
        double psi1 = k * d * sin(thetaAngle);
        double psi2 = k * d * sin(phiAngle);

        // Амплитуда антенны 
        double Af = 1.0; // Амплитуда антенны 
        // Расчет паттерна 
        pattern = Af * (sin(numElementsX * psi1 / 2) / sin(psi1 / 2)) * (sin(numElementsY * psi2 / 2) / sin(psi2 / 2));

        return std::abs(pattern);
    }

    void calculateLOSAngles(const UserTerminal& transmitter, double& phiAOD, double& thetaZOD, double& phiAOA, double& thetaZOA) const {
        double dx = transmitter.x - x;
        double dy = transmitter.y - y;
        double dz = transmitter.z - z;
        double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Углы AOA (угол от приемника к передатчику)
        thetaZOA = acos(dz / distance); // Угловая координата
        phiAOA = atan2(dy, dx); // Азимутальная координата

        // Углы AOD (угол от передатчика к приемнику)
        thetaZOD = acos(-dz / distance); // Угловая координата 
        phiAOD = atan2(-dy, -dx);         // Азимутальная координата 
    }
};

// Функция для вычисления расстояния между двумя пользователями
double calculateDistance(const UserTerminal& a, const UserTerminal& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

class RicianChannel {
public:
    int numPaths;
    double roomWidth;
    double roomLength;
    double roomHeight;

    // delay distribution proportionality factor
    double los_r_tau = 3.0;
    double nlos_r_tau = 3.0;

    std::vector<double> delays_tau;
    

    RicianChannel(int paths, double width, double length, double height)
        : numPaths(paths), roomWidth(width), roomLength(length), roomHeight(height) {}



    // был нужен как черновой вариант ( в модели не используется ) 
    std::complex<double> transmit(const UserTerminal& transmitter, const UserTerminal& receiver, std::complex<double> signal) {
        // Прямое расстояние 
        double directDistance = calculateDistance(transmitter, receiver);
        std::complex<double> directComponent = std::polar(1.0 / directDistance, 0.0); // Прямой компонент

        // Генерация случайных компонентов
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> amplitudeDist(0.0, 1.0);
        std::uniform_real_distribution<> phaseDist(0.0, 2 * M_PI);

        std::complex<double> multipathComponent(0.0, 0.0);
        for (int i = 0; i < numPaths - 1; ++i) {
            double amplitude = amplitudeDist(gen);
            double phase = phaseDist(gen);
            multipathComponent += std::polar(amplitude, phase);
        }

        // Нормализация многопутевого компонента 
        multipathComponent /= std::sqrt(numPaths - 1);

        std::complex<double> receivedSignal = (directComponent + multipathComponent) * signal;

        return receivedSignal;
    }

    // для 3ГГц 
    // Генерация LSP для пользователя для LoS
    LargeScaleParameters generateLargeScaleParametersForLoS() {
        std::random_device rd;
        std::mt19937 gen(rd());

        //Затенение 
        std::lognormal_distribution<> loSshadowFadingDist(0, 3); 
        double loSShadowFading = loSshadowFadingDist(gen);
        
        //К-фактор
        std::normal_distribution<> riceanKDist0(7, 4);
        double riceanK = std::abs(riceanKDist0(gen)); 

        // разброс задержки
        std::normal_distribution<> loSdelaySpreadDist(-7.69802, 0.18);
        double loSDelaySpread = std::abs(loSdelaySpreadDist(gen)); 
        
        //Азимутальный угол Разброс вылета
        std::normal_distribution<> loSAzimuthSpreadDepartureDist(1.6, 0.18);
        double loSAzimuthSpreadDeparture = std::min(std::abs(loSAzimuthSpreadDepartureDist(gen)), 104.0 * M_PI / 180.0);
        
        //Азимутальный угол Разброс прихода
        std::normal_distribution<> loSAzimuthSpreadArrivalDist(1.69037, 0.07224);
        double loSAzimuthSpreadArrival = std::min(std::abs(loSAzimuthSpreadDepartureDist(gen)), 104.0 * M_PI / 180.0);

        //Зенитный угол Разброс прихода     
        double loSZenithSpreadArrival = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.28, 0.239918)), 52.0 * M_PI / 180.0);

        //Зенитный угол Разброс вылета
        double loSZenithSpreadDeparture = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.28, 0.239918)), 52.0 * M_PI / 180.0);

        return LargeScaleParameters(loSShadowFading, riceanK, loSDelaySpread, loSAzimuthSpreadDeparture, loSAzimuthSpreadArrival, loSZenithSpreadDeparture, loSZenithSpreadArrival);
    }

    // Генерация LSP для пользователя для NLoS
    LargeScaleParameters generateLargeScaleParametersForNLoS() {
        std::random_device rd;
        std::mt19937 gen(rd());

        //Затенение
        std::lognormal_distribution<> nLoSshadowFadingDist(0, 8.03); 
        double nLoSShadowFading = nLoSshadowFadingDist(gen);
        //К-фактор = N/A ( заглушка)
        double nLoSRiceanK = 0.0;

        // разброс задержки   
        std::normal_distribution<> nLoSdelaySpreadDist(-7.34158, 0.115206);
        double nLoSDelaySpread = std::abs(nLoSdelaySpreadDist(gen)); // Задержка распространения НЛОС

        //Азимутальный угол Разброс вылета
        std::normal_distribution<> nLoSAzimuthSpreadDepartureDist(1.62, 0.25);
        double nLoSAzimuthSpreadDeparture = std::min(std::abs(nLoSAzimuthSpreadDepartureDist(gen)), 104.0 * M_PI / 180.0);

        //Азимутальный угол Разброс прихода
        std::normal_distribution<> nLoSAzimuthSpreadArrivaDist(1.69037, 0.07224);
        double nLoSAzimuthSpreadArrival = std::min(std::abs(nLoSAzimuthSpreadDepartureDist(gen)), 104.0 * M_PI / 180.0);

        //Зенитный угол Разброс прихода 
        std::normal_distribution<> normalDist4(0, 1);
        double nLoSZenithSpreadArrival = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.28, 0.239918)), 52.0 * M_PI / 180.0);

        //Зенитный угол Разброс вылета
        double nLoSZenithSpreadDeparture = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.28, 0.239918)), 52.0 * M_PI / 180.0);

        return LargeScaleParameters(nLoSShadowFading, nLoSRiceanK, nLoSDelaySpread, nLoSAzimuthSpreadDeparture, nLoSAzimuthSpreadArrival, nLoSZenithSpreadDeparture, nLoSZenithSpreadArrival);
    }




    // Генерация SSP
    std::vector<double> generateClusterDelays(int numClusters, double losdelaySpread, double nlosdelaySpread, double riceanK) {
        std::random_device rd;
        std::mt19937 gen(rd());

        double delay_tau_n;

        for (int n = 0; n < numClusters; ++n) {
            double Xn = std::uniform_real_distribution<>(0.0, 1.0)(gen); // Xn ~ uniform(0,1)

            
            if (n < 15) { 
                delay_tau_n = -1 * los_r_tau * log(Xn) * losdelaySpread ;
                
            }
            else { 
                delay_tau_n = -1 * nlos_r_tau * log(Xn) * nlosdelaySpread ;
            }
            delays_tau.push_back(delay_tau_n);
        }

        // Нормализация задержек
        double minDelay_tau_n = *std::min_element(delays_tau.begin(), delays_tau.end());
        for (auto& delay_tau_n : delays_tau) {
            delay_tau_n -= minDelay_tau_n; // Вычитаем минимальное значение
        }

        // Сортируем задержки
        std::sort(delays_tau.begin(), delays_tau.end());

        /*
        // Если LOS, применяем масштабирование (The scaled delays are not to be used in cluster power generation. ) 
        if (riceanK > 0) { // Если K-фактор положителен
            double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) + (0.0433 *riceanK) + 0.7705;

            delay[0] /= scalingFactor; // Масштабирование задержек

        }
        */

        return delays_tau; // Возвращаем нормализованные задержки
    }


    // Генерация мощностей кластеров
    std::vector<double> generateClusterPowers(const std::vector<double>& clusterDelays, double riceanK, double losdelaySpread, double nlosdelaySpread) {

        riceanK = pow(10, riceanK / 10); // перевод из дб в линейный масштаб

        std::vector<double> clusterPowers(clusterDelays.size());
        double power;
        double maxPower = 0.0;
        double sumclusterPowers = 0.0;
        double losPower = 0.0;

        // Генерация случайных значений для затенения кластеров 
        std::random_device rd;
        std::mt19937 gen(rd());

        std::normal_distribution<> shadowFadingDist(0, 36); // is the per cluster shadowing term in [dB].


        for (size_t n = 0; n < clusterDelays.size(); ++n) {


            double shadowing = shadowFadingDist(gen); // Генерация затенения
            if(n < 15){  
                power = exp( (-1) * clusterDelays[n] * (los_r_tau - 1) / (los_r_tau * losdelaySpread)) * pow(10, (shadowing / 10));
            }
            else{
                power = exp( (-1) * clusterDelays[n] * (nlos_r_tau  - 1) / (nlos_r_tau * nlosdelaySpread)) * pow(10, (shadowing / 10));
            }

            if (n < 15) { // Если LOS, добавляем дополнительный компонент
                losPower = (riceanK / (riceanK + 1));
            }

            clusterPowers[n] = power;

            if (power > maxPower) {
                maxPower = power; // Находим максимальную мощность
            }
        }

        std::cout << "maxPower = " << maxPower - pow(10, 25 / 10) << ": ";

        double threshold = maxPower * 0.03; // Порог для удаления

        auto newEnd = std::remove_if(clusterPowers.begin(), clusterPowers.end(), // убираем значения составляют менее 0,03 от максимума (-25дб) 
            [threshold](double power) { return power < threshold; });
        clusterPowers.erase(newEnd, clusterPowers.end());


        for (auto& n : clusterPowers) sumclusterPowers += n;

        // Нормализация мощностей кластеров
        for (size_t n = 0; n < clusterPowers.size(); ++n) {
            clusterPowers[n] = (1 / (riceanK + 1)) * clusterPowers[n] / sumclusterPowers; // Нормализуем по суммарной мощности
            if (n == 0) { // Если LOS, добавляем дополнительный компонент
                clusterPowers[n] += losPower;
            }
        }

        return clusterPowers; // Возвращаем нормализованные мощности кластеров
    }

    //Scaling factors for AOA, AOD generation
    private:
        double calculateC_phi(double riceanK) {
            // Определяем C_phi в зависимости от K-фактора
            if (riceanK > 0) {
                return 1.211* ((0.0001*riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
            }
            else {
                return 1.273; // NLOS значение
            }
        }
    //Scaling factors for ZOA, ZOD generation
    private:
        double calculateC_theta(double riceanK) {
            // Определяем C_theta в зависимости от K-фактора
            if (riceanK > 0) {
                return 1.1088 * (( 0.0002 * riceanK * riceanK * riceanK) + (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
            }
            else {
                return 1.184; // NLOS значение
            }
        }

    


    // Генерация AOA
    std::vector<double> generateAOA(int numClusters, double AzimuthSpreadArrival, double riceanK) {
        std::vector<double> aoaAngles(numClusters);
        double C_phi = calculateC_phi(riceanK);

        for (int n = 0; n < numClusters; ++n) {
            double Pn = static_cast<double>(n) / (numClusters - 1); // Нормированное значение
            double AOA_n = -C_phi * log(Pn) / (0.4 * AzimuthSpreadArrival); // Применение уравнения (7.5-9)

            // Применение случайного знака
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> signDist(0, 1);
            int sign = (signDist(gen) == 0) ? 1 : -1;

            // Добавление случайного варианта
            AOA_n += sign * generateNormalRandom(0, AzimuthSpreadArrival, gen);
            aoaAngles[n] = AOA_n;
        }

        return aoaAngles;
    }
    /*
    // Генерация AOD (аналогично AOA)
    std::vector<double> generateAOD(int numClusters, double rmsAngleSpread, double riceanK) {
        std::vector<double> aodAngles(numClusters);
        double C_phi = calculateC_phi(riceanK);

        for (int n = 0; n < numClusters; ++n) {
            double Pn = static_cast<double>(n) / (numClusters - 1);
            double AOD_n = -C_phi * log(Pn) / (0.4 * rmsAngleSpread);

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> signDist(0, 1);
            int sign = (signDist(gen) == 0) ? 1 : -1;

            AOD_n += sign * generateNormalRandom(0, rmsAngleSpread, gen);
            aodAngles[n] = AOD_n;
        }

        return aodAngles;
    }

    double generateLaplacianRandom(double mu, double b, std::mt19937& gen) {
        std::laplace_distribution<double> dis(mu, b);
        return dis(gen);
    }

    double calculateC_theta(double riceanK) {
        // Определяем C_theta в зависимости от K-фактора
        if (riceanK > 0) {
            return 0.0002 * pow(riceanK, 3) + 0.027 * pow(riceanK, 2) + 0.15 * riceanK + 0.0339;
        }
        else {
            return 1.3086; // NLOS значение
        }
    }

    std::vector<double> generateZOA(int numClusters, double rmsAngleSpread, double riceanK) {
        std::vector<double> zoaAngles(numClusters);
        double C_theta = calculateC_theta(riceanK);

        for (int n = 0; n < numClusters; ++n) {
            double Pn = static_cast<double>(n) / (numClusters - 1);
            double ZOA_n = -C_theta * log(Pn) / (0.4 * rmsAngleSpread);

            // Случайный знак
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> signDist(0, 1);
            int sign = (signDist(gen) == 0) ? 1 : -1;

            // Добавление случайного варианта
            ZOA_n += sign * generateNormalRandom(0, rmsAngleSpread, gen);
            zoaAngles[n] = ZOA_n;
        }

        return zoaAngles;
    }

    // Генерация ZOD
    std::vector<double> generateZOD(int numClusters, double rmsAngleSpread, double riceanK) {
        std::vector<double> zodAngles(numClusters);
        double C_theta = calculateC_theta(riceanK);

        for (int n = 0; n < numClusters; ++n) {
            double Pn = static_cast<double>(n) / (numClusters - 1);
            double ZOD_n = -C_theta * log(Pn) / (0.4 * rmsAngleSpread);

            // Случайный знак
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> signDist(0, 1);
            int sign = (signDist(gen) == 0) ? 1 : -1;

            // Добавление случайного варианта
            ZOD_n += sign * generateNormalRandom(0, rmsAngleSpread, gen);
            zodAngles[n] = ZOD_n;
        }

        return zodAngles;
    }

    

    



    // Генерация нормального распределения
    double generateNormalRandom(double mean, double stddev, std::mt19937& gen) {
        std::normal_distribution<double> dis(mean, stddev);
        return dis(gen);
    }
    */

};

int main() {
    // Параметры комнаты 
    double roomLength = 120.0; // Длина комнаты
    double roomWidth = 50.0;    // Ширина комнаты 
    double roomHeight = 3.0;    // Высота комнаты 
    double wavelength = 0.3;     // Длина волны 


    // Создание канала Rician
    RicianChannel channel(400, roomWidth, roomLength, roomHeight);

    // Создаем пользователей 
    std::vector<UserTerminal> users;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> xDist(0.0, roomLength);
    std::uniform_real_distribution<> yDist(0.0, roomWidth);
    double userHeight = 1.0; // Высота пользователей 

    for (int i = 0; i < 12; ++i) {
        double bearing = (rand() % 360) * M_PI / 180.0; // Угол поворота 
        double downtilt = (rand() % 90) * M_PI / 180.0; // Угол наклона
        double slant = (rand() % 360) * M_PI / 180.0; // Угол наклона

        users.emplace_back(i, xDist(gen), yDist(gen), userHeight, wavelength, bearing, downtilt, slant);
    }
    // Генерация LSP для каждого пользователя 
    for (auto& user : users) {
        LargeScaleParameters lspForLoS = channel.generateLargeScaleParametersForLoS();
        LargeScaleParameters lspForNLoS = channel.generateLargeScaleParametersForNLoS();
        //std::cout << "User    " << user.id << " LSP: "
            //<< "SF: " << lsp.shadowFading << ", K: " << lsp.riceanK
            //<< ", DS: " << lsp.delaySpread << ", ASA: " << lsp.azimuthSpreadArrival << ", ASD: " << lsp.azimuthSpreadDeparture << ", ZSA: " << lsp.zenithSpreadArrival << ", ZSD: " << lsp.zenithSpreadDeparture << std::endl;

        // Генерация задержек кластеров
        std::vector<double> clusterDelays = channel.generateClusterDelays(34, lspForLoS.delaySpread, lspForNLoS.delaySpread, lspForLoS.riceanK); // Передаем delaySpread из LSP
        //std::cout << "Cluster delays for User " << user.id << ": ";
        for (const auto& delay : clusterDelays) {
            std::cout << delay << "\n ";
        }
        std::cout << std::endl;

        // Генерация мощностей кластеров
        std::vector<double> clusterPowers = channel.generateClusterPowers(clusterDelays, lspForLoS.riceanK, lspForLoS.delaySpread, lspForNLoS.delaySpread);
        std::cout << "Cluster powers for User " << user.id << ": ";
        for (const auto& power : clusterPowers) {
            std::cout << power << " ";
        }
        std::cout << std::endl;
    }
   
    /*
    // Пример использования
    int numClusters = 15; // Количество кластеров
    double rmsAngleSpread = 0.1; // RMS угловой разброс
    double riceanK = 10; // K-фактор в дБ

    std::vector<double> aoaAngles = channel.generateAOA(numClusters, rmsAngleSpread, riceanK);
    std::vector<double> aodAngles = channel.generateAOD(numClusters, rmsAngleSpread, riceanK);

    // Вывод углов AOA и AOD
    std::cout << "AOA Angles: ";
    for (const auto& angle : aoaAngles) {
        std::cout << angle << " ";
    }
    std::cout << std::endl;

    std::cout << "AOD Angles: ";
    for (const auto& angle : aodAngles) {
        std::cout << angle << " ";
    }
    std::cout << std::endl;


    std::vector<double> zoaAngles = channel.generateZOA(numClusters, rmsAngleSpread, riceanK);
    std::vector<double> zodAngles = channel.generateZOD(numClusters, rmsAngleSpread, riceanK);

    // Вывод углов ZOA и ZOD
    std::cout << "ZOA Angles: ";
    for (const auto& angle : zoaAngles) {
        std::cout << angle << " ";
    }
    std::cout << std::endl;

    std::cout << "ZOD Angles: ";
    for (const auto& angle : zodAngles) {
        std::cout << angle << " ";
    }
    std::cout << std::endl;

    */




    // Выбор пар пользователей для передачи
    std::uniform_int_distribution<> userDist(0, users.size() - 1);
    std::set<std::pair<int, int>> selectedPairs;

    while (selectedPairs.size() < 6) {
        int user1 = userDist(gen);
        int user2 = userDist(gen);

        if (user1 != user2) {
            selectedPairs.emplace(std::min(user1, user2), std::max(user1, user2));
        }
    }
















    // Передача сигнала между пользователями 
    std::complex<double> signal(1.0, 0.0); // Сигнал с амплитудой 1 и фазой 0

    for (const auto& pair : selectedPairs) {
        const UserTerminal& transmitter = users[pair.first];
        const UserTerminal& receiver = users[pair.second];

        double phiAOD, thetaZOD, phiAOA, thetaZOA;
        transmitter.calculateLOSAngles(receiver, phiAOD, thetaZOD, phiAOA, thetaZOA);

        std::complex<double> receivedSignal = channel.transmit(transmitter, receiver, signal);


        


        /*
        // Печать информации о передаче 
        double patternTransmitter = transmitter.antennaPattern(thetaZOD, phiAOD);
        double patternReceiver = receiver.antennaPattern(thetaZOA, phiAOA);
        
        std::cout << "Transmission from User " << transmitter.id << " to User " << receiver.id
            << " - Received Signal: " << receivedSignal
            << " | Transmitter Pattern: " << patternTransmitter << " | Receiver Pattern: " << patternReceiver << std::endl;

        std::cout << "LOS AOD: (" << phiAOD << ", " << thetaZOD << "), "
            << "LOS AOA: (" << phiAOA << ", " << thetaZOA << ")\n" << std::endl;

        // Вычисление расстояния между передатчиком и приемником
        double d = calculateDistance(transmitter, receiver);
        std::cout << "Distance between transmitter and receiver is " << d << std::endl;

        double P_LOS;
        if (d <= 5)
            P_LOS = 1;
        else if (d > 5 && d <= 49)
            P_LOS = exp(-(d - 5) / 70.8);
        else
            P_LOS = exp(-(d - 49) / 211.7) * 0.54;

        std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << std::endl << std::endl;*/
    }
    
    return 0;
}
