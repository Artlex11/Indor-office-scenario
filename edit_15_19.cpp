#define _USE_MATH_DEFINES // Для C++
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <set>
#include <algorithm> // Для std::sort
#include <Eigen/Dense>

using namespace Eigen;

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

//___________________________________Класс для STEP_4 Генерация_______________________________________//

class LargeScaleParameters {
public:
    double shadowFading; // Затенение 
    double riceanK; // Фактор Райса 
    double delaySpread; // Задержка распространения 
    double azimuthSpreadDeparture; // Угловой спред выхода
    double azimuthSpreadArrival; // Угловой спред входа 
    double zenithSpreadDeparture; // Зенитный спред выхода 
    double zenithSpreadArrival; // Зенитный спред входа 

    // Конструктор класса
    LargeScaleParameters(double sf, double k, double ds, double asd, double asa, double zsd, double zsa)
        : shadowFading(sf), riceanK(k), delaySpread(ds), azimuthSpreadDeparture(asd),
        azimuthSpreadArrival(asa), zenithSpreadDeparture(zsd), zenithSpreadArrival(zsa) {}
};

//_______________________________________Класс UT________________________________//

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

    // Конструктор класса

    UserTerminal(int id, double x, double y, double z, double lambda, double bearing, double downtilt, double slant)
        : id(id), x(x), y(y), z(z), wavelength(lambda), bearingAngle(bearing), downtiltAngle(downtilt), slantAngle(slant) {
        d = lambda / 2.0; // Расстояние между элементами антенны 
    }

    // Метод для вычисления ДН полей

    double getFieldPatternTheta(double thetaAngle, double phiAngle) const {
        thetaAngle = thetaAngle * 180 / M_PI;
        phiAngle = phiAngle * 180 / M_PI;

        double ksi = +45.0 * 180 / M_PI;
        double f2_Theta = sqrt(std::min(12 * (phiAngle / 65) * (phiAngle / 65), 30.0));
        double f2_Phi = sqrt(std::min(12 * ((thetaAngle - 90) / 65) * ((thetaAngle - 90) / 65), 30.0));

        double cos_Pci = (cos(ksi) * sin(thetaAngle) + sin(ksi) * sin(phiAngle) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));
        double sin_Pci = (sin(ksi) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));

        double f1_Theta = sin_Pci * f2_Theta + cos_Pci * f2_Phi;

        return std::abs(f1_Theta);
    }

    double getFieldPatternPhi(double thetaAngle, double phiAngle) const {
        thetaAngle = thetaAngle * 180 / M_PI;
        phiAngle = phiAngle * 180 / M_PI;

        double ksi = +45.0 * 180 / M_PI;
        double f2_Theta = sqrt(std::min(12 * (phiAngle / 65) * (phiAngle / 65), 30.0));
        double f2_Phi = sqrt(std::min(12 * ((thetaAngle - 90) / 65) * ((thetaAngle - 90) / 65), 30.0));

        double cos_Pci = (cos(ksi) * sin(thetaAngle) + sin(ksi) * sin(phiAngle) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));
        double sin_Pci = (sin(ksi) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));

        double f1_Phi = cos_Pci * f2_Theta - sin_Pci * f2_Phi;

        return std::abs(f1_Phi);
    }




    void calculateLOSAngles(const UserTerminal& transmitter, const UserTerminal& receiver,
        double& losPhiAOD, double& losThetaZOD,
        double& losPhiAOA, double& losThetaZOA) const {

        // Разница координат между передатчиком и приемником
        double dx = receiver.x - transmitter.x;
        double dy = receiver.y - transmitter.y;
        double dz = receiver.z - transmitter.z;
        double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Углы AOA (угол от приемника к передатчику)
        losThetaZOA = acos(dz / distance); // Угловая координата
        losPhiAOA = atan2(dy, dx); // Азимутальная координата

        // Углы AOD (угол от передатчика к приемнику)
        losThetaZOD = acos(-dz / distance); // Угловая координата 
        losPhiAOD = atan2(-dy, -dx);         // Азимутальная координата 
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
    double los_r_tau = 3.6;
    double nlos_r_tau = 3.0;

    std::vector<double> delays_tau;

    // Конструктор класса
    RicianChannel(int paths, double width, double length, double height)
        : numPaths(paths), roomWidth(width), roomLength(length), roomHeight(height) {}




    // для 3ГГц 
    // Генерация LSP(Large Scale Parameters) для пользователя для LoS
    LargeScaleParameters generateLargeScaleParametersForLoS() {
        std::random_device rd;
        std::mt19937 gen(rd());

        //Затенение (SF)
        std::lognormal_distribution<> loSshadowFadingDist(0, 3);
        double loSShadowFading = loSshadowFadingDist(gen);

        //К-фактор (K)
        std::normal_distribution<> riceanKDist0(7, 4);
        double riceanK = std::abs(riceanKDist0(gen));

        // разброс задержки (DS)
        std::normal_distribution<> loSdelaySpreadDist(-7.69802, 0.18);
        double loSDelaySpread = std::abs(loSdelaySpreadDist(gen));

        //Азимутальный угол Разброс вылета (ASD)
        std::normal_distribution<> loSAzimuthSpreadDepartureDist(1.6, 0.18);
        double loSAzimuthSpreadDeparture = std::min(std::abs(loSAzimuthSpreadDepartureDist(gen)), log10(104.0));

        //Азимутальный угол Разброс прихода (ASA)
        std::normal_distribution<> loSAzimuthSpreadArrivalDist(1.69037, 0.07224);
        double loSAzimuthSpreadArrival = std::min(std::abs(loSAzimuthSpreadDepartureDist(gen)), log10(104.0));


        //For frequencies below 6 GHz, use fc = 6 when determining the values of the frequency-dependent ZSD 
        //and ZOD offset parameters in Table 7.5 - 7 and 7.5 - 10

        //Зенитный угол Разброс прихода (ZSA)
        double loSZenithSpreadArrival = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.22027, 0.230196)), log10(52.0));

        // Нет в таблице 7.5-6
        //Зенитный угол Разброс вылета (ZSD)
        double loSZenithSpreadDeparture = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.01951, 0.409863)), log10(52.0));

        return LargeScaleParameters(loSShadowFading, riceanK, loSDelaySpread, loSAzimuthSpreadDeparture, loSAzimuthSpreadArrival, loSZenithSpreadDeparture, loSZenithSpreadArrival);
    }

    // Генерация LSP(Large Scale Parametrs) для пользователя для NLoS
    LargeScaleParameters generateLargeScaleParametersForNLoS() {
        std::random_device rd;
        std::mt19937 gen(rd());

        //Затенение (SF)
        std::lognormal_distribution<> nLoSshadowFadingDist(0, 8.03);
        double nLoSShadowFading = nLoSshadowFadingDist(gen);

        // N/A в таблице
        //К-фактор = N/A (заглушка) (K) 
        double nLoSRiceanK = 0.0;

        // разброс задержки (DS)
        std::normal_distribution<> nLoSdelaySpreadDist(-7.34158, 0.115206);
        double nLoSDelaySpread = std::abs(nLoSdelaySpreadDist(gen)); // Задержка распространения НЛОС

        //Азимутальный угол Разброс вылета (ASD)
        std::normal_distribution<> nLoSAzimuthSpreadDepartureDist(1.62, 0.25);
        double nLoSAzimuthSpreadDeparture = std::min(std::abs(nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));

        //Азимутальный угол Разброс прихода (ASA)
        std::normal_distribution<> nLoSAzimuthSpreadArrivaDist(1.69037, 0.07224);
        double nLoSAzimuthSpreadArrival = std::min(std::abs(nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));


        //Зенитный угол Разброс прихода (ZSA)
        double nLoSZenithSpreadArrival = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.26024, 0.669941)), log10(52.0));

        //Зенитный угол Разброс вылета (ZSD)
        double nLoSZenithSpreadDeparture = std::min(std::abs(laplaceDistribution(generateNormalRandom(0, 1, gen), 1.08, 0.36)), log10(52.0));

        return LargeScaleParameters(nLoSShadowFading, nLoSRiceanK, nLoSDelaySpread, nLoSAzimuthSpreadDeparture, nLoSAzimuthSpreadArrival, nLoSZenithSpreadDeparture, nLoSZenithSpreadArrival);
    }


    //__________________________________________STEP_5____________________________________________________//

    // Генерация SSP(Small Scale Parameters)

    std::vector<double> generateClusterDelays(double losdelaySpread, double nlosdelaySpread, double riceanK)
        {
            std::random_device rn;
            std::mt19937 gen(rn());

            double delay_tau_n;
            delays_tau.clear();

            std::random_device(rd);
            std::default_random_engine re(rd());
            std::uniform_int_distribution<int>dist(0, 1);
            int k = dist(re);
            
            /*delay_tau_n = -1 * los_r_tau * log(Xn) * losdelaySpread;
            delay_tau_n = -1 * nlos_r_tau * log(Xn) * nlosdelaySpread;
            delays_tau.push_back(delay_tau_n);*/

            if (k == 0) {
                for (int n = 0; n < 15; ++n)
                {
                    double Xn = std::uniform_real_distribution<>(0.0, 1.0)(gen); // Xn ~ uniform(0,1)
                    delay_tau_n = -1 * los_r_tau * log(Xn) * losdelaySpread;
                    delays_tau.push_back(delay_tau_n);
                }
            }
            else {
                for (int n = 0; n < 19; ++n)
                {
                    double Xn = std::uniform_real_distribution<>(0.0, 1.0)(gen); // Xn ~ uniform(0,1)
                    delay_tau_n = -1 * nlos_r_tau * log(Xn) * nlosdelaySpread;
                    delays_tau.push_back(delay_tau_n);
                }
                
            }
            
    

            // Нормализация задержек
            double minDelay_tau_n = *std::min_element(delays_tau.begin(), delays_tau.end());
            for (auto& delay_tau_n : delays_tau) {
                delay_tau_n -= minDelay_tau_n; // Вычитаем минимальное значение
            }

            // Сортируем задержки
            std::sort(delays_tau.begin(), delays_tau.end());

            // дополнительный параметр для LOS лучей, вычисляется через K-фактор, дальше не учитывается
            /*
            // Если LOS, применяем масштабирование (The scaled delays are not to be used in cluster power generation. )
            if (riceanK > 0) { // Если K-фактор положителен
                double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) + (0.0433 *riceanK) + 0.7705;

                delay[0] /= scalingFactor; // Масштабирование задержек

            }
            */

            return delays_tau; // Возвращаем нормализованные задержки
    }

    //__________________________________________STEP_6____________________________________________________//

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

        for (size_t n = 0; n < clusterDelays.size(); ++n) {
           
            if (n == 15) {
                std::normal_distribution<> shadowFadingDist(0, 36); // is the per cluster shadowing term in [dB] LOS.
                double shadowing = shadowFadingDist(gen); // Генерация затенения
                power = exp((-1) * clusterDelays[n] * (los_r_tau - 1) / (los_r_tau * losdelaySpread)) * pow(10, (shadowing / 10));
            }
            else {
                std::normal_distribution<> shadowFadingDist(0, 9); // is the per cluster shadowing term in [dB] NLOS.
                double shadowing = shadowFadingDist(gen); // Генерация затенения
                power = exp((-1) * clusterDelays[n] * (nlos_r_tau - 1) / (nlos_r_tau * nlosdelaySpread)) * pow(10, (shadowing / 10));
            }

            if (n == 0) { // Если LOS, добавляем дополнительный компонент
                losPower = (riceanK / (riceanK + 1));
            }

            clusterPowers[n] = power;
        }

        // Нормализация мощностей кластеров
        maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
        for (auto& n : clusterPowers) sumclusterPowers += n;

        for (size_t n = 0; n < clusterPowers.size(); ++n) {
            clusterPowers[n] = (1 / (riceanK + 1)) * clusterPowers[n] / sumclusterPowers; // Нормализуем по суммарной мощности
            if (n == 0) { // Если LOS, добавляем дополнительный компонент
                clusterPowers[n] += losPower;
            }
        }

        double threshold = maxPower / sumclusterPowers * 0.003; // Порог для удаления
        for (size_t n = 0; n < clusterPowers.size(); ++n) {
            if (clusterPowers[n] < threshold) { clusterPowers[n] = 0; }
        }

        return clusterPowers; // Возвращаем нормализованные мощности кластеров
    }

    //______________________________________________STEP_9________________________________________________//

    MatrixXd generateXPR(int numClusters) {

        std::random_device rd;
        std::mt19937 gen(rd());

        MatrixXd XPR(numClusters, 20);
        for (int n = 0; n < numClusters; ++n) {
            for (int m = 0; m < 20; ++m) {
                std::normal_distribution<> X_n_m_Dist(1.62, 0.25);
                XPR(n, m) = pow(10, (X_n_m_Dist(gen)) / 10);
            }
        }
        return XPR;
    }

    //______________________________________________STEP_10________________________________________________//

    // Метод для генерации случайных начальных фаз
    MatrixXd generateInitialRandomPhases(int numClusters, int numRay) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> phaseDist(-M_PI, M_PI); // Распределение фаз от -π до π

        MatrixXd initialRandomPhases(numClusters, numPaths * 4);

        for (int n = 0; n < numClusters; ++n) {
            for (int m = 0; m < numRay; ++m) {
                // Генерация случайных фаз для каждой комбинации поляризации
                initialRandomPhases(n, m * 4) = phaseDist(gen); // Theta - Theta
                initialRandomPhases(n, m * 4 + 1) = phaseDist(gen); // Theta - Phi
                initialRandomPhases(n, m * 4 + 2) = phaseDist(gen); // Phi - Theta
                initialRandomPhases(n, m * 4 + 3) = phaseDist(gen); // Phi - Phi
            }
        }
        return initialRandomPhases;
    }




    //______________________________________________STEP_7________________________________________________//

    // Cluster ASA , ASD , ZSA in [deg] for los and nlos
    double los_C_ASA = 8.0;
    double los_C_ASD = 5.0;
    double los_C_ZSA = 9.0;
    double los_C_ZSD = 0.375;

    double nlos_C_ASA = 11.0;
    double nlos_C_ASD = 5.0;
    double nlos_C_ZSA = 9.0;
    double nlos_C_ZSD = 0.375;

    std::vector<double> am = { 0.0447, -0.0447 ,0.1413 ,-0.1413,0.2492 ,-0.2492 ,0.3715 ,-0.3715 ,0.5129 ,-0.5129 ,
        0.6797 ,-0.6797, 0.8844 , -0.8844,1.1481,-1.1481, 1.5195 , -1.5195, 2.1551 , -2.1551 };



    //AOA
    MatrixXd generatePhiAOA(const std::vector<double>& clusterPowers, int numClusters, double losAzimuthSpreadArrival, double nlosAzimuthSpreadArrival, double riceanK, double losPhiAOA) {

        std::random_device rd;
        std::mt19937 gen(rd());


        double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
        std::vector<double> Phi_n_AOA(clusterPowers.size());
        MatrixXd Phi_n_m_AOA(clusterPowers.size(), 20);


        double X1;
        double Y1;
        double Phi_1_AOA;


        for (int n = 0; n < numClusters; ++n) {

            double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

            if (n < 15) {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                Phi_n_AOA[n] = (2 * (losAzimuthSpreadArrival / 1.4) * pow(-log(clusterPowers[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
                double Yn = generateNormalRandom(0, (losAzimuthSpreadArrival / 7) * (losAzimuthSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Phi_1_AOA = Phi_n_AOA[n]; }

                Phi_n_AOA[n] = Phi_n_AOA[n] * Xn + Yn + losPhiAOA - Phi_1_AOA * X1 - Y1;


                for (int m = 0; m < 20; ++m) {
                    Phi_n_m_AOA(n, m) = Phi_n_AOA[n] + los_C_ASA * am[m];
                }
            }
            else {
                double C_phi = 1.273;
                Phi_n_AOA[n] = (2 * (nlosAzimuthSpreadArrival / 1.4) * pow(-log(clusterPowers[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
                double Yn = generateNormalRandom(0, (nlosAzimuthSpreadArrival / 7) * (nlosAzimuthSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Phi_1_AOA = Phi_n_AOA[n]; }

                Phi_n_AOA[n] = Phi_n_AOA[n] * Xn + Yn + losPhiAOA;

                for (int m = 0; m < 20; ++m) {
                    Phi_n_m_AOA(n, m) = Phi_n_AOA[n] + nlos_C_ASA * am[m];
                }
            }
        }

        return Phi_n_m_AOA;
    }
    //AOD
    MatrixXd generatePhiAOD(const std::vector<double>& clusterPowers, int numClusters, double losAzimuthSpreadDeparture, double nlosAzimuthSpreadDeparture, double riceanK, double losPhiAOD) {

        std::random_device rd;
        std::mt19937 gen(rd());


        double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
        std::vector<double> Phi_n_AOD(clusterPowers.size());
        MatrixXd Phi_n_m_AOD(clusterPowers.size(), 20);


        double X1;
        double Y1;
        double Phi_1_AOD;


        for (int n = 0; n < numClusters; ++n) {

            double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

            if (n < 15) {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                Phi_n_AOD[n] = (2 * (losAzimuthSpreadDeparture / 1.4) * pow(-log(clusterPowers[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (losAzimuthSpreadDeparture / 7) * (losAzimuthSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Phi_1_AOD = Phi_n_AOD[n]; }

                Phi_n_AOD[n] = Phi_n_AOD[n] * Xn + Yn + losPhiAOD - Phi_1_AOD * X1 - Y1;


                for (int m = 0; m < 20; ++m) {
                    Phi_n_m_AOD(n, m) = Phi_n_AOD[n] + los_C_ASD * am[m];
                }
            }
            else {
                double C_phi = 1.273;
                Phi_n_AOD[n] = (2 * (nlosAzimuthSpreadDeparture / 1.4) * pow(-log(clusterPowers[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (nlosAzimuthSpreadDeparture / 7) * (nlosAzimuthSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Phi_1_AOD = Phi_n_AOD[n]; }

                Phi_n_AOD[n] = Phi_n_AOD[n] * Xn + Yn + losPhiAOD;


                for (int m = 0; m < 20; ++m) {
                    Phi_n_m_AOD(n, m) = Phi_n_AOD[n] + nlos_C_ASD * am[m];
                }
            }
        }

        return Phi_n_m_AOD;
    }

    //ZOA
    MatrixXd generateThetaZOA(const std::vector<double>& clusterPowers, int numClusters, double losZenithSpreadArrival, double nlosZenithSpreadArrival, double riceanK, double losThetaZOA) {

        std::random_device rd;
        std::mt19937 gen(rd());

        double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
        std::vector<double> Theta_n_ZOA(clusterPowers.size());
        MatrixXd Theta_n_m_ZOA(clusterPowers.size(), 20);


        double X1;
        double Y1;
        double Theta_1_ZOA;

        for (int n = 0; n < numClusters; ++n) {


            double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

            if (n < 15) {
                double C_theta = 1.184 * ((0.0002 * riceanK * riceanK * riceanK) + (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
                Theta_n_ZOA[n] = ((-1) * losZenithSpreadArrival * (clusterPowers[n] / maxPower)) / C_theta; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (losZenithSpreadArrival / 7) * (losZenithSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Theta_1_ZOA = Theta_n_ZOA[n]; }

                Theta_n_ZOA[n] = Theta_n_ZOA[n] * Xn + Yn + losThetaZOA - Theta_1_ZOA * X1 - Y1;


                for (int m = 0; m < 20; ++m) {
                    Theta_n_m_ZOA(n, m) = Theta_n_ZOA[n] + los_C_ZSA * am[m];
                }

            }
            else {
                double C_theta = 1.184;
                Theta_n_ZOA[n] = ((-1) * nlosZenithSpreadArrival * (clusterPowers[n] / maxPower)) / C_theta; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (nlosZenithSpreadArrival / 7) * (nlosZenithSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Theta_1_ZOA = Theta_n_ZOA[n]; }

                Theta_n_ZOA[n] = Theta_n_ZOA[n] * Xn + Yn + losThetaZOA;

                for (int m = 0; m < 20; ++m) {
                    Theta_n_m_ZOA(n, m) = Theta_n_ZOA[n] + nlos_C_ZSA * am[m];
                }
            }
        }
        return Theta_n_m_ZOA;
    }

    // ZOD
    MatrixXd generateThetaZOD(const std::vector<double>& clusterPowers, int numClusters, double losZenithSpreadDeparture, double nlosZenithSpreadDeparture, double riceanK, double losThetaZOD) {

        std::random_device rd;
        std::mt19937 gen(rd());

        double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
        std::vector<double> Theta_n_ZOD(clusterPowers.size());
        MatrixXd Theta_n_m_ZOD(clusterPowers.size(), 20);

        double X1;
        double Y1;
        double Theta_1_ZOD;

        for (int n = 0; n < numClusters; ++n) {

            double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

            if (n < 15) {
                double C_theta = 1.184 * ((0.0002 * riceanK * riceanK * riceanK) + (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
                Theta_n_ZOD[n] = ((-1) * losZenithSpreadDeparture * (clusterPowers[n] / maxPower)) / C_theta; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (losZenithSpreadDeparture / 7) * (losZenithSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Theta_1_ZOD = Theta_n_ZOD[n]; }


                Theta_n_ZOD[n] = Theta_n_ZOD[n] * Xn + Yn + losThetaZOD - Theta_1_ZOD * X1 - Y1;

                for (int m = 0; m < 20; ++m) {
                    Theta_n_m_ZOD(n, m) = Theta_n_ZOD[n] + los_C_ZSD * am[m];
                }
            }
            else {
                double C_theta = 1.184;
                Theta_n_ZOD[n] = ((-1) * nlosZenithSpreadDeparture * (clusterPowers[n] / maxPower)) / C_theta; // Применение уравнения (7.5-9)
                double Yn = generateNormalRandom(0, (nlosZenithSpreadDeparture / 7) * (nlosZenithSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

                if (n == 0) { X1 = Xn; Y1 = Yn; Theta_1_ZOD = Theta_n_ZOD[n]; }

                Theta_n_ZOD[n] = Theta_n_ZOD[n] * Xn + Yn + losThetaZOD;

                for (int m = 0; m < 20; ++m) {
                    Theta_n_m_ZOD(n, m) = Theta_n_ZOD[n] + nlos_C_ZSD * am[m];
                }
            }
        }
        return Theta_n_m_ZOD;
    }
};


//_______________________________________Основная_Программа___________________________________________//

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

    for (int i = 1; i <= 12; ++i) {
        double bearing = (rand() % 360) * M_PI / 180.0; // Угол поворота 
        double downtilt = (rand() % 90) * M_PI / 180.0; // Угол наклона
        double slant = (rand() % 360) * M_PI / 180.0; // Угол наклона

        users.emplace_back(i, xDist(gen), yDist(gen), userHeight, wavelength, bearing, downtilt, slant);
    }


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



        // Расчёт лос угов и соовтетствующих значений ДН полей приёмника и передатчика
        double losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA;
        transmitter.calculateLOSAngles(transmitter, receiver, losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA);

        std::cout << "LOS AOD ZOD: (" << losPhiAOD << ", " << losThetaZOD << "), "
            << "LOS AOA ZOA: (" << losPhiAOA << ", " << losThetaZOA << ")\n" << std::endl;

        double F_tx_Theta = transmitter.getFieldPatternTheta(losThetaZOD, losPhiAOD);
        double F_tx_Pfi = transmitter.getFieldPatternPhi(losThetaZOD, losPhiAOD);
        std::cout << "F_tx_Theta(theta,pfi) : " << F_tx_Theta << "\n ";
        std::cout << "F_tx_Pfi(theta,pfi) : " << F_tx_Pfi << "\n ";

        double F_rx_Theta = receiver.getFieldPatternTheta(losThetaZOD, losPhiAOD);
        double F_rx_Pfi = receiver.getFieldPatternPhi(losThetaZOA, losPhiAOA);
        std::cout << "F_rx_Theta(theta,pfi) : " << F_rx_Theta << "\n ";
        std::cout << "F_rx_Pfi(theta,pfi) : " << F_rx_Pfi << "\n ";

        LargeScaleParameters lspForLoS = channel.generateLargeScaleParametersForLoS();
        LargeScaleParameters lspForNLoS = channel.generateLargeScaleParametersForNLoS();

        // Вывод LSP параметров для LOS
        std::cout << "LSP for LOS for User " << transmitter.id << "and User " << receiver.id << " LSP: "
            << "SF: " << lspForLoS.shadowFading << ", K: " << lspForLoS.riceanK
            << ", DS: " << lspForLoS.delaySpread << ", ASA: " << lspForLoS.azimuthSpreadArrival << ", ASD: " << lspForLoS.azimuthSpreadDeparture << ", ZSA: " << lspForLoS.zenithSpreadArrival << ", ZSD: " << lspForLoS.zenithSpreadDeparture << std::endl << std::endl;

        // Вывод LSP параметров для NLOS
        std::cout << "LSP for NLOS for User " << transmitter.id << "and User " << receiver.id << " LSP: "
            << "SF: " << lspForNLoS.shadowFading << ", K: " << lspForNLoS.riceanK
            << ", DS: " << lspForNLoS.delaySpread << ", ASA: " << lspForNLoS.azimuthSpreadArrival << ", ASD: " << lspForNLoS.azimuthSpreadDeparture << ", ZSA: " << lspForNLoS.zenithSpreadArrival << ", ZSD: " << lspForNLoS.zenithSpreadDeparture << std::endl;

        //______________STEP_5_______________//
        // Генерация задержек кластеров
        std::vector<double> clusterDelays = channel.generateClusterDelays(lspForLoS.delaySpread, lspForNLoS.delaySpread, lspForLoS.riceanK); // Передаем delaySpread из LSP
        std::cout << "Cluster delays for User " << transmitter.id << "and User " << receiver.id << ":" << std::endl;
        int i = 1;
        for (const auto& delay : clusterDelays) {
            std::cout << i << "-delay: " << delay << "\n ";
            i++;
        }
        std::cout << std::endl;

        //_____________STEP_6_______________//
        // Генерация мощностей кластеров
        std::vector<double> clusterPowers = channel.generateClusterPowers(clusterDelays, lspForLoS.riceanK, lspForLoS.delaySpread, lspForNLoS.delaySpread);
        std::cout << "Cluster powers for User " << transmitter.id << "and User " << receiver.id << ": ";
        for (const auto& power : clusterPowers) {
            std::cout << power << " ";
        }
        std::cout << std::endl << std::endl;

        //_____________STEP_7_______________//

        //AOA
        MatrixXd PhiAOA = channel.generatePhiAOA(clusterPowers, clusterDelays.size(), lspForLoS.azimuthSpreadArrival, lspForNLoS.azimuthSpreadArrival, lspForLoS.riceanK, losPhiAOA);
        std::cout << "PhiAOA for UT receiver " << receiver.id << ": \n";
        // Вывод значений матрицы построчно
        for (int i = 0; i < PhiAOA.rows(); ++i) {
            for (int j = 0; j < PhiAOA.cols(); ++j) {
                std::cout << PhiAOA(i, j) << " "; // Вывод элемента
            }
            std::cout << std::endl; // Переход на новую строку после вывода всех элементов строки
        }
        std::cout << std::endl;

        //AOD
        MatrixXd PhiAOD = channel.generatePhiAOD(clusterPowers, clusterDelays.size(), lspForLoS.azimuthSpreadDeparture, lspForNLoS.azimuthSpreadDeparture, lspForLoS.riceanK, losPhiAOD);
        std::cout << "PhiAOD for UT transmitter " << transmitter.id << ": \n";
        // Вывод значений матрицы построчно
        for (int i = 0; i < PhiAOD.rows(); ++i) {
            for (int j = 0; j < PhiAOD.cols(); ++j) {
                std::cout << PhiAOD(i, j) << " "; // Вывод элемента
            }
            std::cout << std::endl; // Переход на новую строку после вывода всех элементов строки
        }
        std::cout << std::endl;

        //ZOA
        MatrixXd ThetaZOA = channel.generateThetaZOA(clusterPowers, clusterDelays.size(), lspForLoS.zenithSpreadArrival, lspForNLoS.azimuthSpreadArrival, lspForLoS.riceanK, losThetaZOA);
        std::cout << "ThetaZOA for UT receiver " << receiver.id << ": \n";
        // Вывод значений матрицы построчно
        for (int i = 0; i < ThetaZOA.rows(); ++i) {
            for (int j = 0; j < ThetaZOA.cols(); ++j) {
                std::cout << ThetaZOA(i, j) << " "; // Вывод элемента
            }
            std::cout << std::endl; // Переход на новую строку после вывода всех элементов строки
        }
        std::cout << std::endl;

        //ZOD
        MatrixXd ThetaZOD = channel.generateThetaZOD(clusterPowers, clusterDelays.size(), lspForLoS.zenithSpreadDeparture, lspForNLoS.azimuthSpreadDeparture, lspForLoS.riceanK, losThetaZOD);
        std::cout << "ThetaZOD for UT transmitter " << transmitter.id << ": \n";
        // Вывод значений матрицы построчно
        for (int i = 0; i < ThetaZOD.rows(); ++i) {
            for (int j = 0; j < ThetaZOD.cols(); ++j) {
                std::cout << ThetaZOD(i, j) << " "; // Вывод элемента
            }
            std::cout << std::endl; // Переход на новую строку после вывода всех элементов строки
        }
        std::cout << std::endl;

        //_____________STEP_9_______________//

        //XPR
        MatrixXd XPR = channel.generateXPR(clusterDelays.size());
        std::cout << "Generate the cross polarization power ratios  K_n_m: \n";
        for (int i = 0; i < XPR.rows(); ++i) {
            for (int j = 0; j < XPR.cols(); ++j) {
                std::cout << XPR(i, j) << " "; // Вывод элемента
            }
            std::cout << std::endl; // Переход на новую строку после вывода всех элементов строки
        }
        std::cout << std::endl;

        //_____________STEP_10_______________//
        /*MatrixXd initialPhases = channel.generateInitialRandomPhases(clusterDelays.size(), 20);
        
        std::cout << "Initial random phases for each ray in each cluster:" << std::endl;
        for (int n = 0; n < 34; ++n) {
            std::cout << "Cluster " << n + 1 << ": \n";
            for (int m = 0; m < 20; ++m) {
                std::cout << "Ray " << m + 1 << " (Theta - Theta, Theta - Phi, Phi - Theta, Phi - Phi): | "
                    << initialPhases(n, m * 4) << ", "
                    << initialPhases(n, m * 4 + 1) << ", "
                    << initialPhases(n, m * 4 + 2) << ", "
                    << initialPhases(n, m * 4 + 3) << " |\n ";
            }
            std::cout << std::endl;
        }*/
        
        //________________________________STEP_2___________________________________________________________//

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

        std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << std::endl << std::endl;
       
        //__________________________________STEP_3____________________________________________________________//
        // Вычисление Pathloss_LOS для каждой пары
        double path_loss_LOS;
        double nu = 3; // частота в ГГц
        path_loss_LOS = 32.4 + 17.3 * log10(d) + 20 * log10(nu);

        std::cout << "PathLoss in LOS between users " << transmitter.id << " and " << receiver.id << " = " << path_loss_LOS << std::endl << std::endl;
    }


    return 0;
}


