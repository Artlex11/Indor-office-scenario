#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <set>
#include <algorithm> 
#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <string>
#include <fstream>
#include <filesystem>
#include <iomanip>

using namespace Eigen;

struct GaussGenerator
{
    GaussGenerator(double mean, double stddev, std::uint32_t seed)
        : engine_(seed), distribution_(mean, stddev) {}

    GaussGenerator(double mean, double stddev)
        : distribution_(mean, stddev)
    {
        using namespace std;
        seed_seq seeds{
            (uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count(),
            (uint64_t)chrono::system_clock::now().time_since_epoch().count(),
            (uint64_t)hash<thread::id>{}(this_thread::get_id()),
        };
        engine_.seed(seeds);
    }

    double operator()() { return distribution_(engine_); }

    std::mt19937 engine_;
    std::normal_distribution<double> distribution_;
};

template<typename T>
struct UniformGenerator
{
    UniformGenerator(T min, T max, std::uint32_t seed)
        : engine_(seed), distribution_(min, max) {}

    UniformGenerator(T min, T max)
        : distribution_(min, max)
    {
        using namespace std;
        seed_seq seeds{
            (uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count(),
            (uint64_t)chrono::system_clock::now().time_since_epoch().count(),
            (uint64_t)hash<thread::id>{}(this_thread::get_id()),
        };
        engine_.seed(seeds);
    }

    T operator()() { return distribution_(engine_); }

private:
    std::mt19937 engine_;
    std::conditional_t<std::is_integral<T>::value, std::uniform_int_distribution<T>, std::uniform_real_distribution<T>> distribution_;
};


template<typename T>
class UniformGeneratorFromVector {
public:
    // Конструктор, принимающий вектор значений
    UniformGeneratorFromVector(const std::vector<T>& values)
        : values_(values) {
        // Генерация семени
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine_.seed(seeds);
        // Инициализация распределения для выбора индекса
        distribution_ = std::uniform_int_distribution<size_t>(0, values_.size() - 1);
    }

    // Оператор вызова для получения случайного элемента
    T operator()() {
        return values_[distribution_(engine_)];
    }

private:
    std::mt19937 engine_; // Генератор случайных чисел
    std::vector<T> values_; // Вектор значений
    std::uniform_int_distribution<size_t> distribution_; // Распределение для выбора индекса
};


std::vector<int> indicesToDelete;

//_________________________Класс_для_STEP_4_Генерация_LSP______________________________________//
class LargeScaleParameters {
public:
    double shadowFading; // Затенение 
    double riceanK; // Фактор Райса 
    double delaySpread; // Задержка распространения 
    double azimuthSpreadDeparture; // Угловой спред выхода
    double azimuthSpreadArrival; // Угловой спред входа 
    double zenithSpreadDeparture; // Зенитный спред выхода 
    double zenithSpreadArrival; // Зенитный спред входа 
    bool los;
    double fc;

    // Конструктор для лос
    LargeScaleParameters(bool losvalue, double fcvalue) : los(losvalue), fc(fcvalue)
    {
        double StandardDeviationSF;
        double StandardDeviationK;
        double StandardDeviationDS;
        double StandardDeviationASD;
        double StandardDeviationASA;
        double StandardDeviationZSA;
        double StandardDeviationZSD;

        double MeanSF;
        double MeanK;
        double MeanDS;
        double MeanASD;
        double MeanASA;
        double MeanZSA;
        double MeanZSD;

        if (los)
        {
            VectorXd value(7);
            VectorXd means(7);
            MatrixXd C(7, 7);

            //Если LOS :
            //Затемнение SF
            StandardDeviationSF = 3;
            MeanSF = 0;
            //К-фактор (K)
            StandardDeviationK = 4;
            MeanK = 7;
            // разброс задержки (DS)
            StandardDeviationDS = 0.18;
            MeanDS = (-0.01 * log10(1 + fc) - 7.692);
            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.18;
            MeanASD = 1.6;
            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.119);
            MeanASA = (-0.19 * log10(1 + fc) + 1.781);
            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.04 * log10(1 + fc) + 0.264);
            MeanZSA = (-0.26 * log10(1 + fc) + 1.44);
            //Зенитный угол Разброс вылета (ZSD)
            if (fc < 6.0) {
                StandardDeviationZSD = (0.13 * log10(1 + 6) + 0.30);
                MeanZSD = (-1.43 * log10(1 + 6) + 2.228);
            }
            else {
                StandardDeviationZSD = (0.13 * log10(1 + fc) + 0.30);
                MeanZSD = (-1.43 * log10(1 + fc) + 2.228);
            }

            value << shadowFading, riceanK, delaySpread, azimuthSpreadDeparture, azimuthSpreadArrival, zenithSpreadDeparture, zenithSpreadArrival;
            means << MeanSF, MeanK, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

            //SF K DS ASD ASA ZSD ZSA
            C << StandardDeviationSF * StandardDeviationSF, 0.5 * StandardDeviationSF * StandardDeviationK, -0.8 * StandardDeviationSF * StandardDeviationDS, -0.4 * StandardDeviationSF * StandardDeviationASD, -0.5 * StandardDeviationSF * StandardDeviationASA, 0.2 * StandardDeviationSF * StandardDeviationZSD, 0.3 * StandardDeviationSF * StandardDeviationZSA,
                0.5 * StandardDeviationK * StandardDeviationSF, StandardDeviationK* StandardDeviationK, -0.5 * StandardDeviationK * StandardDeviationDS, 0.0 * StandardDeviationK * StandardDeviationASD, 0.0 * StandardDeviationK * StandardDeviationASA, 0.0 * StandardDeviationK * StandardDeviationZSD, 0.1 * StandardDeviationK * StandardDeviationZSA,
                -0.8 * StandardDeviationDS * StandardDeviationSF, -0.5 * StandardDeviationDS * StandardDeviationK, StandardDeviationDS* StandardDeviationDS, 0.6 * StandardDeviationDS * StandardDeviationASD, 0.8 * StandardDeviationDS * StandardDeviationASA, 0.1 * StandardDeviationDS * StandardDeviationZSD, 0.2 * StandardDeviationDS * StandardDeviationZSA,
                -0.4 * StandardDeviationASD * StandardDeviationSF, 0.0 * StandardDeviationASD * StandardDeviationK, 0.6 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.4 * StandardDeviationASD * StandardDeviationASA, 0.5 * StandardDeviationASD * StandardDeviationZSD, 0.0 * StandardDeviationASD * StandardDeviationZSA,
                -0.5 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationK, 0.8 * StandardDeviationASA * StandardDeviationDS, 0.4 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, 0.0 * StandardDeviationASA * StandardDeviationZSD, 0.5 * StandardDeviationASA * StandardDeviationZSA,
                0.2 * StandardDeviationZSD * StandardDeviationSF, 0.0 * StandardDeviationZSD * StandardDeviationK, 0.1 * StandardDeviationZSD * StandardDeviationDS, 0.5 * StandardDeviationZSD * StandardDeviationASD, 0.0 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.0 * StandardDeviationZSD * StandardDeviationZSA,
                0.3 * StandardDeviationZSA * StandardDeviationSF, 0.1 * StandardDeviationZSA * StandardDeviationK, 0.2 * StandardDeviationZSA * StandardDeviationDS, 0.0 * StandardDeviationZSA * StandardDeviationASD, 0.5 * StandardDeviationZSA * StandardDeviationASA, 0.0 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

            MatrixXd L;
            L = C.llt().matrixL();

            GaussGenerator rand(0.0, 1.0);
            for (int i = 0; i < 7; ++i) {
                value(i) = rand();
            }

            VectorXd value_new = L * value + means;

            shadowFading = value_new(0);
            riceanK = value_new(1);
            delaySpread = value_new(2);
            azimuthSpreadDeparture = std::min((value_new(3)), log10(104.0));
            azimuthSpreadArrival = std::min((value_new(4)), log10(104.0));
            zenithSpreadDeparture = std::min((value_new(5)), log10(52.0));
            zenithSpreadArrival = std::min((value_new(6)), log10(52.0));

        }
        else
        {
            VectorXd value(6);
            VectorXd means(6);
            MatrixXd C(6, 6);

            //Если NLOS:
            //Затенение (SF)
            StandardDeviationSF = 8.03;
            MeanSF = 0;
            // разброс задержки (DS)
            StandardDeviationDS = (0.10 * log10(1 + fc) + 0.055);
            MeanDS = (-0.28 * log10(1 + fc) - 7.173);
            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.25;
            MeanASD = 1.62;
            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.059);
            MeanASA = (-0.11 * log10(1 + fc) + 1.863);
            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.09 * log10(1 + fc) + 0.746);
            MeanZSA = (-0.15 * log10(1 + fc) + 1.387);
            //Зенитный угол Разброс вылета (ZSD)
            StandardDeviationZSD = 0.36;
            MeanZSD = 1.08;

            value << shadowFading, delaySpread, azimuthSpreadDeparture, azimuthSpreadArrival, zenithSpreadDeparture, zenithSpreadArrival;
            means << MeanSF, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

            //SF  DS ASD ASA ZSD ZSA
            C << StandardDeviationSF * StandardDeviationSF, -0.5 * StandardDeviationSF * StandardDeviationDS, 0.0 * StandardDeviationSF * StandardDeviationASD, -0.4 * StandardDeviationSF * StandardDeviationASA, 0.0 * StandardDeviationSF * StandardDeviationZSD, 0.0 * StandardDeviationSF * StandardDeviationZSA,
                -0.5 * StandardDeviationDS * StandardDeviationSF, StandardDeviationDS* StandardDeviationDS, 0.4 * StandardDeviationDS * StandardDeviationASD, 0.0 * StandardDeviationDS * StandardDeviationASA, -0.27 * StandardDeviationDS * StandardDeviationZSD, -0.06 * StandardDeviationDS * StandardDeviationZSA,
                0.0 * StandardDeviationASD * StandardDeviationSF, 0.4 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.0 * StandardDeviationASD * StandardDeviationASA, 0.35 * StandardDeviationASD * StandardDeviationZSD, 0.23 * StandardDeviationASD * StandardDeviationZSA,
                -0.4 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationDS, 0.0 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, -0.08 * StandardDeviationASA * StandardDeviationZSD, 0.43 * StandardDeviationASA * StandardDeviationZSA,
                0.0 * StandardDeviationZSD * StandardDeviationSF, -0.27 * StandardDeviationZSD * StandardDeviationDS, 0.35 * StandardDeviationZSD * StandardDeviationASD, -0.08 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.42 * StandardDeviationZSD * StandardDeviationZSA,
                0.0 * StandardDeviationZSA * StandardDeviationSF, -0.06 * StandardDeviationZSA * StandardDeviationDS, 0.23 * StandardDeviationZSA * StandardDeviationASD, 0.43 * StandardDeviationZSA * StandardDeviationASA, 0.42 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

            MatrixXd L;
            L = C.llt().matrixL();
            GaussGenerator rand(0.0, 1.0);
            for (int i = 0; i < 6; ++i) {
                value(i) = rand();
            }

            // Преобразование случайных величин с помощью матрицы Лапласа
            VectorXd value_new = L * value + means;

            shadowFading = value_new(0);
            delaySpread = value_new(1);
            azimuthSpreadDeparture = std::min((value_new(2)), log10(104.0));
            azimuthSpreadArrival = std::min((value_new(3)), log10(104.0));
            zenithSpreadDeparture = std::min((value_new(4)), log10(52.0));
            zenithSpreadArrival = std::min((value_new(5)), log10(52.0));
        }
    }

    void showParameters() {
        if (los) {
            // Вывод LSP параметров для NLOS
            std::cout << "SF [dB] : " << shadowFading << ",\nK [dB] : " << riceanK << ",\nDS [log10(DS/1s)] : " << delaySpread << ",\nASD [log10(ASA/ 1* degree] : " << azimuthSpreadDeparture << ",\nASA [log10(ASA/ 1* degree] :  " << azimuthSpreadArrival << ",\nZSD [log10(ZSD/ 1* degree] : " << zenithSpreadDeparture << ",\nZSA [log10(ZSA/ 1* degree] : " << zenithSpreadArrival << std::endl << std::endl;
        }
        else {
            // Вывод LSP параметров для NLOS
            std::cout << "SF [dB] : " << shadowFading << ",\nDS [log10(DS/1s)] : " << delaySpread << ",\nASD [log10(ASD/ 1* degree] : " << azimuthSpreadDeparture << ",\nASA [log10(ASA/ 1* degree] : " << azimuthSpreadArrival << ",\nZSD [log10(ZSD/ 1* degree] : " << zenithSpreadDeparture << ",\nZSA [log10(ZSA/ 1* degree] : " << zenithSpreadArrival << std::endl << std::endl;
        }
    }

};

//__________________________________________________Класс_UT______________________________________________//
class UserTerminal {
public:
    int id;
    double x, y, z; // Координаты пользователя 
    double wavelength = 0.1; // Длина волны 
    int numElementsX = 4; // 4 элемента антенны по X
    int numElementsY = 2; // 2 элемента антенны по Y 
    double bearingAngle; // Угол поворота 
    double downtiltAngle; // Угол наклона 
    double slantAngle; // Угол наклона

    // Конструктор класса(пустой)

    UserTerminal(int id, double x, double y, double z, double bearing, double downtilt, double slant)
        : id(id), x(x), y(y), z(z), bearingAngle(bearing), downtiltAngle(downtilt), slantAngle(slant)
    {}

    // Методы для вычисления ДН полей //Vector2d - вектор столбцы с двумя координатами
    Vector2d FieldPattern(double theta, double phi, double ksi) const
    {

        double alpha = bearingAngle;
        double beta = downtiltAngle;
        double gamma = slantAngle;

        //std::cout << alpha << ", " << beta << ", " << gamma;
        // Inverse rotation matrix
        MatrixXd inv_R(3, 3);
        inv_R <<
            cos(alpha) * cos(beta), sin(alpha)* cos(beta), -sin(beta),
            cos(alpha)* sin(beta)* sin(gamma) - sin(alpha) * cos(gamma), sin(alpha)* sin(beta)* sin(gamma) + cos(alpha) * cos(gamma), cos(beta)* sin(gamma),
            cos(alpha)* sin(beta)* cos(gamma) + sin(alpha) * sin(gamma), sin(alpha)* sin(beta)* cos(gamma) - cos(alpha) * sin(gamma), cos(beta)* cos(gamma);
        //std::cout << inv_R;

        VectorXd xyz(3);
        xyz << sin(theta * M_PI / 180) * cos(phi * M_PI / 180), sin(theta * M_PI / 180)* sin(phi * M_PI / 180), cos(theta * M_PI / 180); // r = 1 
        Vector3d v(0, 0, 1);
        Vector3cd u(1.0, std::complex<double>(0.0, 1.0), 0.0);


        double theta_LSC = acos(v.transpose() * inv_R * xyz);
        std::complex<double> temp = u.transpose() * inv_R * xyz;
        //std::cout << temp1;
        double phi_LSC = atan2(temp.imag(), temp.real());



        //while (phi < -180) {
        //    phi += 360; // Добавляем 360, пока значение не станет >= -180
        //}
        //while (phi > 180) {
        //    phi -= 360; // Вычитаем 360, пока значение не станет <= 180
        //}

        //// Нормализация theta в диапазоне [0, 180]
        //while (theta < 0) {
        //    theta += 180; // Добавляем 180, пока значение не станет >= 0
        //}
        //while (theta > 180) {
        //    theta -= 180; // Вычитаем 180, пока значение не станет <= 180
        //}

        Eigen::Vector2d fieldPattern;
        phi_LSC = phi_LSC * 180 / M_PI;
        theta_LSC = theta_LSC * 180 / M_PI;
        //std::cout << theta_LSC << ", " << phi_LSC << std::endl;
        ksi = ksi * M_PI / 180;



        // Theta
        // -min(12((Theta-90)/Theta3D)^2, SLA_V) Theta3D = 65 deg SLA_V = 30 dB
        double AverticalPowerPattern = (-1) * std::min(12 * ((theta_LSC - 90) / 65) * ((theta_LSC - 90) / 65), 30.0);

        // Phi
        //-min(12(Phi/Phi3D)^2, A_max) Phi_3D = 65 deg, A_max = 30 dB
        double AhorizontalPowerPattern = (-1) * std::min(12 * (phi_LSC / 65) * (phi_LSC / 65), 30.0);



        double A3D_PowerPattern = -std::min((-1) * (AverticalPowerPattern + AhorizontalPowerPattern), 30.0); // знак "-" в начале формулы

        A3D_PowerPattern = pow(10, A3D_PowerPattern / 20);


        double f1_Theta = sqrt(A3D_PowerPattern) * cos(ksi);
        double f1_Phi = sqrt(A3D_PowerPattern) * sin(ksi);

        fieldPattern << f1_Theta, f1_Phi;

        return fieldPattern;
    }

    // в fieldPattern - 1 theta, 2 phi

    void calculateLOSAngles(const UserTerminal& transmitter, const UserTerminal& receiver,
        double& losPhiAOD, double& losThetaZOD,
        double& losPhiAOA, double& losThetaZOA) const
    {
        double dx = receiver.x - transmitter.x;
        double dy = receiver.y - transmitter.y;
        double dz = receiver.z - transmitter.z;

        double distance3D = std::sqrt(dx * dx + dy * dy + dz * dz);
        //double distance2D = std::sqrt(dx * dx + dy * dy );

        losThetaZOA = acos(-dz / distance3D); 
        losPhiAOA = atan2(-dy, -dx); 

        losThetaZOD = acos(dz / distance3D); 
        losPhiAOD = atan2(dy, dx);        
    }

    /*

                                ^Z
                                |
                / (0)   / (1)   |   / (2)  / (3)
                                |
                / (4)   / (5)   |   / (6)  / (7)
                                X- - - - - - - - - - - - - - - - - ->y
                \ (0)   \ (1)       \ (2)  \ (3)

                \ (4)   \ (5)       \ (6)  \ (7)


    */

    MatrixXd generateAntennaElements() const
    {
        double alpha = bearingAngle;
        double beta = downtiltAngle;
        double gamma = slantAngle;

        MatrixXd R(3, 3);
        R <<
            cos(alpha) * cos(beta), cos(alpha)* sin(beta)* sin(gamma) - sin(alpha) * cos(gamma), cos(alpha)* sin(beta)* cos(gamma) + sin(alpha) * sin(gamma),
            sin(alpha)* cos(beta), sin(alpha)* sin(beta)* sin(gamma) + cos(alpha) * cos(gamma), sin(alpha)* sin(beta)* cos(gamma) - cos(alpha) * sin(gamma),
            -sin(beta), cos(beta)* sin(gamma), cos(beta)* cos(gamma);



        MatrixXd locationMatrixAntennaElements(16, 3);

        locationMatrixAntennaElements.row(0) << 0, -3 * wavelength / 4, wavelength / 4;
        locationMatrixAntennaElements.row(8) = locationMatrixAntennaElements.row(0);
        locationMatrixAntennaElements.row(1) << 0, -wavelength / 4, wavelength / 4;
        locationMatrixAntennaElements.row(9) = locationMatrixAntennaElements.row(1);
        locationMatrixAntennaElements.row(2) << 0, wavelength / 4, wavelength / 4;
        locationMatrixAntennaElements.row(10) = locationMatrixAntennaElements.row(2);
        locationMatrixAntennaElements.row(3) << 0, 3 * wavelength / 4, wavelength / 4;
        locationMatrixAntennaElements.row(11) = locationMatrixAntennaElements.row(3);
        locationMatrixAntennaElements.row(4) << 0, -3 * wavelength / 4, -wavelength / 4;
        locationMatrixAntennaElements.row(12) = locationMatrixAntennaElements.row(4);
        locationMatrixAntennaElements.row(5) << 0, -wavelength / 4, -wavelength / 4;
        locationMatrixAntennaElements.row(13) = locationMatrixAntennaElements.row(5);
        locationMatrixAntennaElements.row(6) << 0, wavelength / 4, -wavelength / 4;
        locationMatrixAntennaElements.row(14) = locationMatrixAntennaElements.row(6);
        locationMatrixAntennaElements.row(7) << 0, 3 * wavelength / 4, -wavelength / 4;
        locationMatrixAntennaElements.row(15) = locationMatrixAntennaElements.row(7);

        locationMatrixAntennaElements = (R * locationMatrixAntennaElements.transpose()).transpose();

        return locationMatrixAntennaElements;
    }
};

//___________________Функция_для_вычисления_расстояния_между_двумя_пользователями____________________//
double calculateDistance(const UserTerminal& a, const UserTerminal& b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}


//______________________________________________STEP_5_____________________________________________________//
//________________________________Генерация_SSP(Small_Scale_Parameters)____________________________________//
//______________________________delay_distribution_proportionality_factor__________________________________//
//_______________________________Задаю_один_раз,_нужны_в_STEP_5_и_STEP_6___________________________________//
double los_r_tau = 3.6;
double nlos_r_tau = 3.0;

std::vector<double> generateClusterDelays(bool los, double delaySpread, double riceanK)
{
    delaySpread = pow(10, delaySpread);
    UniformGenerator<double> XnDist(0.0, 1.0);

    std::vector<double> delays_tau;
    double delay_tau_n;
    //delays_tau.clear();

    if (los) {
        for (int n = 0; n < 15; ++n)
        {
            double Xn = XnDist(); // Xn ~ uniform(0,1)
            delay_tau_n = -1 * los_r_tau * log(Xn) * delaySpread;
            delays_tau.push_back(delay_tau_n);
        }
    }
    else {
        for (int n = 0; n < 19; ++n)
        {
            double Xn = XnDist(); // Xn ~ uniform(0,1)
            delay_tau_n = -1 * nlos_r_tau * log(Xn) * delaySpread;
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
    
    // Если LOS, применяем масштабирование (The scaled delays are not to be used in cluster power generation. )
    //if (los) { // Если K-фактор положителен
    //    double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) + (0.0433 *riceanK) + 0.7705;
    //    for (int n = 0; n < delays_tau.size(); ++n)
    //    {
    //       delays_tau[n] /= scalingFactor; // Масштабирование задержек
    //    }
    //}
    
    return delays_tau; // Возвращаем нормализованные задержки
};

//_______________________________________________STEP_6____________________________________________________//
//__________________________________Генерация мощностей кластеров__________________________________________//
std::vector<double> generateClusterPowers(bool los, const std::vector<double>& clusterDelays, double delaySpread)
{
    delaySpread = pow(10, delaySpread);

    std::vector<double> clusterPowers(clusterDelays.size());
    double power = 0.0;
    double maxPower = 0.0;
    double sumclusterPowers = 0.0;

    // Генерация случайных значений для затенения кластеров 
    

    for (size_t n = 0; n < clusterDelays.size(); ++n) {

        if (los) {
            GaussGenerator shadowFadingDist(0.0, 6.0); 
            double shadowing = shadowFadingDist(); 
            power = exp(((-1) * clusterDelays[n] * (los_r_tau - 1))  / (los_r_tau * delaySpread)) * pow(10, (-0.1 * shadowing));
        }
        else {
            GaussGenerator shadowFadingDist(0.0, 3.0); 
            double shadowing = shadowFadingDist(); 
            power = exp(((-1) * clusterDelays[n] * (nlos_r_tau - 1))  / (nlos_r_tau * delaySpread)) * pow(10, (-0.1 * shadowing));
        }
        clusterPowers[n] = power;
    }

    // Нормализация мощностей кластеров
    for (auto& n : clusterPowers) sumclusterPowers += n;

    for (size_t n = 0; n < clusterPowers.size(); ++n) {
        clusterPowers[n] = clusterPowers[n] / sumclusterPowers;  
    }
    
    return clusterPowers; // Возвращаем нормализованные мощности кластеров
};

//________________________________________________STEP_7___________________________________________________//
//____________________________Cluster_ASA_,_ASD_,_ZSA_in_[deg]_for_los_and_nlos____________________________//
double los_C_ASA = 8.0;
double los_C_ASD = 5.0;
double los_C_ZSA = 9.0;
double los_C_ZSD = 0.375 * pow(10, -1.43 * log(1 + 30) + 2.228) ;

double nlos_C_ASA = 11.0;
double nlos_C_ASD = 5.0;
double nlos_C_ZSA = 9.0;
double nlos_C_ZSD = 0.375 * pow(10, 1.08);

std::vector<double> los_C = { los_C_ASA , los_C_ASD , los_C_ZSA , los_C_ZSD };
std::vector<double> nlos_C = { nlos_C_ASA , nlos_C_ASD , nlos_C_ZSA , nlos_C_ZSD };




std::vector<double> am{ 0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492,
    0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
    0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
    2.1551, -2.1551 };


//_________________________________________________AOA/AOD____________________________________________________//
Eigen::MatrixXd generateAOAorAOD_n_m(bool los, const std::vector<double>& clusterPowers, double ASAorASD, double riceanK, double AOAorAOD, int AOA_0_or_AOD_1) {

    AOAorAOD = AOAorAOD * 180 / M_PI;
    ASAorASD = pow(10, ASAorASD);
    GaussGenerator YnDist(0.0, (ASAorASD / 7.0) );
    UniformGeneratorFromVector <int> XnDist(std::vector<int>{ -1, 1 });

    double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());
    double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) - (0.028 * riceanK) + 1.1035);
    double X1 ;
    double Y1 ;
    double AOAorAOD_n1;

    if (los)
    {
        VectorXd AOAorAOD_n(clusterPowers.size());
        MatrixXd AOAorAOD_n_m(clusterPowers.size() + 1, 20);
        //double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) - (0.028 * riceanK) + 1.1035);
        //AOAorAOD_n_m(0, 0) = AOAorAOD;
        for (int n = 0; n < clusterPowers.size(); n++)
        {
            double Xn = XnDist();
            double Yn = YnDist();

            AOAorAOD_n(n) = 2.0 * (ASAorASD / 1.4) * sqrt(-log(clusterPowers[n] / maxPower)) / C_phi;

            if (n == 0)
            {
                X1 = Xn;
                Y1 = Yn; 
                AOAorAOD_n1 = AOAorAOD_n(n);
                AOAorAOD_n_m(n, 0) = AOAorAOD;

            }
            
            AOAorAOD_n(n) = AOAorAOD_n(n) * Xn + Yn + AOAorAOD - AOAorAOD_n1 * X1 - Y1;

            for (int m = 0; m < 20; ++m)
            {
                AOAorAOD_n_m(n + 1, m) = AOAorAOD_n(n) + los_C[AOA_0_or_AOD_1] * am[m];
                //while (AOAorAOD_n_m(n + 1, m) < -180) {
                //    AOAorAOD_n_m(n + 1, m) += 360; // Добавляем 360, пока значение не станет >= -180
                //}
                //while (AOAorAOD_n_m(n + 1, m) > 180) {
                //    AOAorAOD_n_m(n + 1, m) -= 360; // Вычитаем 360, пока значение не станет <= 180
                //}
            }
        }
        return AOAorAOD_n_m;
    }
    else // случай NLOS
    {
        VectorXd AOAorAOD_n(clusterPowers.size());
        MatrixXd AOAorAOD_n_m(clusterPowers.size(), 20);
        C_phi = 1.273;

        for (int n = 0; n < clusterPowers.size(); n++)
        {
            double Xn = XnDist(); 
            double Yn = YnDist();

            AOAorAOD_n(n) = 2.0 * (ASAorASD / 1.4) * sqrt(-log(clusterPowers[n] / maxPower)) / C_phi;

            AOAorAOD_n(n) = AOAorAOD_n(n) * Xn + Yn + AOAorAOD;

            for (int m = 0; m < 20; ++m)
            {
                
                AOAorAOD_n_m(n, m) = AOAorAOD_n(n) + nlos_C[AOA_0_or_AOD_1] * am[m];
                //while (AOAorAOD_n_m(n , m) < -180) {
                //    AOAorAOD_n_m(n , m) += 360; // Добавляем 360, пока значение не станет >= -180
                //}
                //while (AOAorAOD_n_m(n , m) > 180) {
                //    AOAorAOD_n_m(n , m) -= 360; // Вычитаем 360, пока значение не станет <= 180
                //}
            }
        }
        return AOAorAOD_n_m;
    }
}

//_____________________________________________ZOA/ZOD____________________________________________________//
Eigen::MatrixXd generateZOAorZOD_n_m(bool los, const std::vector<double>& clusterPowers, double ZSAorZSD, double riceanK, double ZOAorZOD, int ZOA_2_or_ZOD_3) {

    ZOAorZOD = ZOAorZOD * 180 / M_PI;
    ZSAorZSD = pow(10, ZSAorZSD);

    GaussGenerator YnDist(0.0, (ZSAorZSD / 7.0) );

    UniformGeneratorFromVector <int> XnDist(std::vector<int>{ -1, 1 });

    double maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());


    double C_theta = 1.184 * ((0.0002 * riceanK * riceanK * riceanK) - (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
    double X1 ;
    double Y1 ;
    double ZOAorZOD_n1;
    double meanOffsetZOD = 0.0;

    if (los) // случай LOS
    {
        VectorXd ZOAorZOD_n(clusterPowers.size());
        MatrixXd ZOAorZOD_n_m(clusterPowers.size() + 1, 20);
        //double C_theta = 1.1088 * ((0.0002 * riceanK * riceanK * riceanK) - (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
        //ZOAorZOD_n_m(0, 0) = ZOAorZOD;
        for (int n = 0; n < clusterPowers.size(); n++)
        {         
            double Xn = XnDist();
                
            double Yn = YnDist();

            ZOAorZOD_n(n) = -1 * ZSAorZSD * log(clusterPowers[n] / maxPower) / C_theta;

            if (n == 0)
            {
                X1 = Xn;
                Y1 = Yn;
                ZOAorZOD_n1 = ZOAorZOD_n(n);
                ZOAorZOD_n_m(n, 0) = ZOAorZOD;
            }
            
            ZOAorZOD_n(n) = ZOAorZOD_n(n) * Xn + Yn + ZOAorZOD - ZOAorZOD_n1 * X1 - Y1 + meanOffsetZOD;

            for (int m = 0; m < 20; ++m)
            {
                ZOAorZOD_n_m(n + 1, m) = ZOAorZOD_n(n) + los_C[ZOA_2_or_ZOD_3] * am[m];

                if ( ZOAorZOD_n_m(n + 1, m) >= 180 && ZOAorZOD_n_m(n + 1, m) <= 360)  {
                    ZOAorZOD_n_m(n + 1, m) = 360 - ZOAorZOD_n_m(n + 1, m); 
                }
            }
            
        }
        return ZOAorZOD_n_m;
    }
    else // случай NLOS
    {
        VectorXd ZOAorZOD_n(clusterPowers.size());
        MatrixXd ZOAorZOD_n_m(clusterPowers.size(), 20);
        C_theta = 1.184;
        for (int n = 0; n < clusterPowers.size(); n++)
        {
            double Xn = XnDist();
            double Yn = YnDist();

            ZOAorZOD_n(n) = -1 * ZSAorZSD * log(clusterPowers[n] / maxPower) / C_theta;

            ZOAorZOD_n(n) = ZOAorZOD_n(n) * Xn + Yn + ZOAorZOD + meanOffsetZOD;

            for (int m = 0; m < 20; ++m)
            {
                ZOAorZOD_n_m(n, m) = ZOAorZOD_n(n) + nlos_C[ZOA_2_or_ZOD_3] * am[m];

                if (ZOAorZOD_n_m(n , m) >= 180 && ZOAorZOD_n_m(n , m) <= 360) {
                    ZOAorZOD_n_m(n , m) = 360 - ZOAorZOD_n_m(n , m);
                }
            }
        }
        return ZOAorZOD_n_m;
    }
}

//_________________________________________STEP_8_____________________________________________________//
// Функция для перемешивания лучей (путей) в 4 углах для каждого кластера
void randomCouplingRays(Eigen::MatrixXd& matrix1, Eigen::MatrixXd& matrix2,
    Eigen::MatrixXd& matrix3, Eigen::MatrixXd& matrix4,
    bool los) {
    std::vector<std::mt19937> gens = {
        std::mt19937(std::random_device()()),
        std::mt19937(std::random_device()()),
        std::mt19937(std::random_device()()),
        std::mt19937(std::random_device()())
    };
    std::vector<Eigen::MatrixXd*> matrices = { &matrix1, &matrix2, &matrix3, &matrix4 };

    //индексы для суб-кластеров
    std::vector<int> category1 = { 0, 1, 2, 3, 4, 5, 6, 7, 18, 19 }; // 1й суб-кластер
    std::vector<int> category2 = { 8, 9, 10, 11, 16, 17 }; // 2й суб-кластер
    std::vector<int> category3 = { 12, 13, 14, 15 }; // 3й суб-кластер

    for (size_t i = 0; i < matrices.size(); ++i) {
        for (int j = 0; j < matrices[i]->rows(); ++j) {
            std::vector<double> row(20); // Размер вектора на 20 элементов
            for (int k = 0; k < 20 && k < matrices[i]->cols(); ++k) {
                row[k] = (*matrices[i])(j, k);
            }
            if (los) {
                // Если los == true
                if (j == 0) {
                    // Первая строка не изменяется
                    continue;
                }
                else if (j == 1 || j == 2) {
                    // Перемешиваем элементы из каждого суб-кластера
                    std::vector<double> selectedElements1, selectedElements2, selectedElements3;
                    for (int index : category1) {
                        selectedElements1.push_back(row[index]);
                    }
                    for (int index : category2) {
                        selectedElements2.push_back(row[index]);
                    }
                    for (int index : category3) {
                        selectedElements3.push_back(row[index]);
                    }
                    // Перемешиваем каждый субкластер
                    std::shuffle(selectedElements1.begin(), selectedElements1.end(), gens[i]);
                    std::shuffle(selectedElements2.begin(), selectedElements2.end(), gens[i]);
                    std::shuffle(selectedElements3.begin(), selectedElements3.end(), gens[i]);
                    // Записываем перемешанные элементы обратно в строку
                    for (size_t k = 0; k < selectedElements1.size(); ++k) {
                        row[category1[k]] = selectedElements1[k];
                    }
                    for (size_t k = 0; k < selectedElements2.size(); ++k) {
                        row[category2[k]] = selectedElements2[k];
                    }
                    for (size_t k = 0; k < selectedElements3.size(); ++k) {
                        row[category3[k]] = selectedElements3[k];
                    }
                }
            }
            else {
                // Если los == false
                if (j == 0 || j == 1) {
                    std::vector<double> selectedElements1, selectedElements2, selectedElements3;
                    for (int index : category1) {
                        selectedElements1.push_back(row[index]);
                    }
                    for (int index : category2) {
                        selectedElements2.push_back(row[index]);
                    }
                    for (int index : category3) {
                        selectedElements3.push_back(row[index]);
                    }
                    // Перемешиваем каждую подкатегорию
                    std::shuffle(selectedElements1.begin(), selectedElements1.end(), gens[i]);
                    std::shuffle(selectedElements2.begin(), selectedElements2.end(), gens[i]);
                    std::shuffle(selectedElements3.begin(), selectedElements3.end(), gens[i]);
                    // Записываем перемешанные элементы обратно в строку
                    for (size_t k = 0; k < selectedElements1.size(); ++k) {
                        row[category1[k]] = selectedElements1[k];
                    }
                    for (size_t k = 0; k < selectedElements2.size(); ++k) {
                        row[category2[k]] = selectedElements2[k];
                    }
                    for (size_t k = 0; k < selectedElements3.size(); ++k) {
                        row[category3[k]] = selectedElements3[k];
                    }
                }
            }
            // Перемешиваем оставшиеся строки (начиная с 3-й строки, если los == true, и с 4-й, если los == false)
            if (j >= (los ? 3 : 2)) {
                std::shuffle(row.begin(), row.begin() + 20, gens[i]); // Перемешиваем только первые 20 элементов
            }
            for (int k = 0; k < 20 && k < matrices[i]->cols(); ++k) {
                (*matrices[i])(j, k) = row[k];
            }
        }
    }
}

//_________________________________________A_1_A_2________________________________________________________//
std::vector<double> calculateAngularSpreadandMeanAngles(bool los, const std::vector<double>& clasterPowers, MatrixXd& AOD, MatrixXd& AOA, MatrixXd& ZOD, MatrixXd& ZOA) {
    std::vector<double> ASandMeanAnglesforAOD_AOA_ZOD_ZOA;

    AOD = AOD * M_PI / 180;
    AOA = AOA * M_PI / 180;
    ZOD = ZOD * M_PI / 180;
    ZOA = ZOA * M_PI / 180;

    std::complex<double> weighted_sumAOA(0.0, 0.0);
    std::complex<double> weighted_sumAOD(0.0, 0.0);
    std::complex<double> weighted_sumZOA(0.0, 0.0);
    std::complex<double> weighted_sumZOD(0.0, 0.0);
    double weighted_sumPowers = 0.0;

    // Вычисление взвешенной суммы комплексных экспонент
    for (int n = 0; n < clasterPowers.size(); ++n) {
        for (int m = 0; m < 20; ++m) {
            if (los) {
                weighted_sumAOD += (clasterPowers[n] / 20) * std::complex<double>(cos(AOD(n + 1, m)), sin(AOD(n + 1, m)));
                weighted_sumAOA += (clasterPowers[n] / 20) * std::complex<double>(cos(AOA(n + 1, m)), sin(AOA(n + 1, m)));
                weighted_sumZOD += (clasterPowers[n] / 20) * std::complex<double>(cos(ZOD(n + 1, m)), sin(ZOD(n + 1, m)));
                weighted_sumZOA += (clasterPowers[n] / 20) * std::complex<double>(cos(ZOA(n + 1, m)), sin(ZOA(n + 1, m)));
                weighted_sumPowers += (clasterPowers[n] / 20);
            }
            else if(!los) {
                weighted_sumAOD += (clasterPowers[n] / 20) * std::complex<double>(cos(AOD(n , m)), sin(AOD(n , m)));
                weighted_sumAOA += (clasterPowers[n] / 20) * std::complex<double>(cos(AOA(n , m)), sin(AOA(n , m)));
                weighted_sumZOD += (clasterPowers[n] / 20) * std::complex<double>(cos(ZOD(n , m)), sin(ZOD(n , m)));
                weighted_sumZOA += (clasterPowers[n] / 20) * std::complex<double>(cos(ZOA(n , m)), sin(ZOA(n , m)));
                weighted_sumPowers += (clasterPowers[n] / 20);
            }


        }
    }

    // Вычисление углового рассеяния
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(sqrt(-2.0 * std::log(std::abs(weighted_sumAOD / weighted_sumPowers))));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(sqrt(-2.0 * std::log(std::abs(weighted_sumAOA / weighted_sumPowers))));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(sqrt(-2.0 * std::log(std::abs(weighted_sumZOD / weighted_sumPowers))));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(sqrt(-2.0 * std::log(std::abs(weighted_sumZOA / weighted_sumPowers))));


    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(atan2(weighted_sumAOD.imag(), weighted_sumAOD.real()));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(atan2(weighted_sumAOA.imag(), weighted_sumAOA.real()));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(atan2(weighted_sumZOD.imag(), weighted_sumZOD.real()));
    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(atan2(weighted_sumZOA.imag(), weighted_sumZOA.real()));

    return ASandMeanAnglesforAOD_AOA_ZOD_ZOA;
}

//_________________________________________________STEP_2__________________________________________//
bool generateLOSorNLOS(double distance , int Open0orMixed1) {
    UniformGenerator<double> pDist(0.0, 1.0); // равномерное распределение вероятности LOS (подброс монетка)
    double P_LOS;
    bool los;
    double p; // случайно сгенерированная вероятность (порог LOS)

    if (Open0orMixed1 == 0) {
        if (distance <= 5.0) {
            P_LOS = 1;
            los = true;
            std::cout << "Link - LOS" << std::endl;
        }
        else if (distance > 5 && distance <= 49) {
            P_LOS = exp(-(distance - 5) / 70.8);
            p = pDist();
            // Сравнение с порогом, если P_LOS > p => LOS, иначе NLOS
            if (P_LOS > p) {
                los = true;
                std::cout << "Link - LOS" << std::endl;
            }
            else {
                los = false;
                std::cout << "Link - NLOS" << std::endl;
            }
        }
        else {
            P_LOS = exp(-(distance - 49) / 211.7) * 0.54;
            p = pDist();
            // Сравнение с порогом, если P_LOS > p => LOS, иначе NLOS
            if (P_LOS > p) {
                los = true;
                std::cout << "Link - LOS" << std::endl;
            }
            else {
                los = false;
                std::cout << "Link - NLOS" << std::endl;
            }
        }
    }

    else if (Open0orMixed1 == 1) {
        if (distance <= 1.2) {
            P_LOS = 1;
            los = true;
            std::cout << "Link - LOS" << std::endl;
        }
        else if (distance > 1.2 && distance <= 6.5) {
            P_LOS = exp(-(distance - 1.2) / 4.7);
            p = pDist();
            // Сравнение с порогом, если P_LOS > p => LOS, иначе NLOS
            if (P_LOS > p) {
                los = true;
                std::cout << "Link - LOS" << std::endl;
            }
            else {
                los = false;
                std::cout << "Link - NLOS" << std::endl;
            }
        }
        else {
            P_LOS = exp(-(distance - 6.5) / 32.6) * 0.32;
            p = pDist();
            // Сравнение с порогом, если P_LOS > p => LOS, иначе NLOS
            if (P_LOS > p) {
                los = true;
                std::cout << "Link - LOS" << std::endl;
            }
            else {
                los = false;
                std::cout << "Link - NLOS" << std::endl;
            }
        }
    }
    return los;
}

//_________________________________________________STEP_3__________________________________________//
double calculatePathLoss(bool los, double distance, double nu) {
    double path_loss;
    if (los)
    {
        path_loss = 32.4 + 17.3 * log10(distance) + 20 * log10(nu);
    }
    else
    {
        path_loss = 38.3 * log10(distance) + 17.30 + 24.9 * log10(nu);
    }
    return path_loss;
}


int main()
{   
    std::ofstream AS_AOA_CDF, AS_AOD_CDF, AS_ZOA_CDF, AS_ZOD_CDF;
    std::ofstream Mean_AOA_CDF, Mean_AOD_CDF, Mean_ZOA_CDF, Mean_ZOD_CDF;

    AS_AOA_CDF << std::fixed << std::setprecision(5);
    AS_AOD_CDF << std::fixed << std::setprecision(5);
    AS_ZOA_CDF << std::fixed << std::setprecision(5);
    AS_ZOD_CDF << std::fixed << std::setprecision(5);

    Mean_AOA_CDF << std::fixed << std::setprecision(5);
    Mean_AOD_CDF << std::fixed << std::setprecision(5);
    Mean_ZOA_CDF << std::fixed << std::setprecision(5);
    Mean_ZOD_CDF << std::fixed << std::setprecision(5);
    

    AS_AOA_CDF.open("AS_AOA_CDF.txt");
    AS_AOD_CDF.open("AS_AOD_CDF.txt");
    AS_ZOA_CDF.open("AS_ZOA_CDF.txt");
    AS_ZOD_CDF.open("AS_ZOD_CDF.txt");

    Mean_AOA_CDF.open("Mean_AOA_CDF.txt");
    Mean_AOD_CDF.open("Mean_AOD_CDF.txt");
    Mean_ZOA_CDF.open("Mean_ZOA_CDF.txt");
    Mean_ZOD_CDF.open("Mean_ZOD_CDF.txt");

    std::ofstream SF_LOS_CDF, K_LOS_CDF, DS_LOS_CDF, ASA_LOS_CDF, ASD_LOS_CDF, ZSA_LOS_CDF, ZSD_LOS_CDF;
    std::ofstream SF_NLOS_CDF, DS_NLOS_CDF, ASA_NLOS_CDF, ASD_NLOS_CDF, ZSA_NLOS_CDF, ZSD_NLOS_CDF;

    SF_LOS_CDF << std::fixed << std::setprecision(10);
    K_LOS_CDF << std::fixed << std::setprecision(10);
    DS_LOS_CDF << std::fixed << std::setprecision(10);
    ASA_LOS_CDF << std::fixed << std::setprecision(10);
    ASD_LOS_CDF << std::fixed << std::setprecision(10);
    ZSA_LOS_CDF << std::fixed << std::setprecision(10);
    ZSD_LOS_CDF << std::fixed << std::setprecision(10);


    SF_NLOS_CDF << std::fixed << std::setprecision(10);
    DS_NLOS_CDF << std::fixed << std::setprecision(10);
    ASA_NLOS_CDF << std::fixed << std::setprecision(10);
    ASD_NLOS_CDF << std::fixed << std::setprecision(10);
    ZSA_NLOS_CDF << std::fixed << std::setprecision(10);
    ZSD_NLOS_CDF << std::fixed << std::setprecision(10);

    
    SF_LOS_CDF.open("SF_LOS_CDF.txt");
    K_LOS_CDF.open("K_LOS_CDF.txt");
    DS_LOS_CDF.open("DS_LOS_CDF.txt");
    ASA_LOS_CDF.open("ASA_LOS_CDF.txt");
    ASD_LOS_CDF.open("ASD_LOS_CDF.txt");
    ZSA_LOS_CDF.open("ZSA_LOS_CDF.txt");
    ZSD_LOS_CDF.open("ZSD_LOS_CDF.txt");

    SF_NLOS_CDF.open("SF_NLOS_CDF.txt");
    DS_NLOS_CDF.open("DS_NLOS_CDF.txt");
    ASA_NLOS_CDF.open("ASA_NLOS_CDF.txt");
    ASD_NLOS_CDF.open("ASD_NLOS_CDF.txt");
    ZSA_NLOS_CDF.open("ZSA_NLOS_CDF.txt");
    ZSD_NLOS_CDF.open("ZSD_NLOS_CDF.txt");


   
                
    // Параметры комнаты 
    double roomLength = 120.0; // Длина комнаты
    double roomWidth = 80.0;    // Ширина комнаты 
    double roomHeight = 3.0;    // Высота комнаты 

    double nu = 30.0; // частота в ГГц
    double wavelength = 30000000 / nu;     // Длина волны 

    
    for (int rooms = 0; rooms <= 1000; ++rooms) {

        // Создаем пользователей 
        std::vector<UserTerminal> users;
        UniformGenerator<double> yDist(0.0, roomLength - 0.0);
        UniformGenerator<double> xDist(0.0, roomWidth - 0.0);
        double userHeight = 1.0; // Высота пользователей 

        for (int i = 1; i <= 12; ++i) {
            double bearing = (rand() % 360) * M_PI / 180.0;
            double downtilt = (rand() % 180) * M_PI / 180.0; // Угол наклона
            double slant = (rand() % 360) * M_PI / 180.0; // Угол наклона

            UserTerminal newUT(i, 0, 0, userHeight, bearing, downtilt, slant);
            bool isValidPosition = false;

            while (!isValidPosition) {
                newUT.x = xDist();
                newUT.y = yDist();

                // Проверяем расстояние до всех существующих пользователей
                isValidPosition = true; // Предполагаем, что позиция валидна
                for (const auto& user : users) {
                    double distanseUTnew_UT = calculateDistance(newUT, user);
                    if (distanseUTnew_UT < 1.5 ) {
                        isValidPosition = false; // Позиция невалидна
                        break; // Выходим из цикла
                    }
                }
                users.push_back(newUT);
            }
        }

        // Выбор пар пользователей для передачи

        UniformGenerator<int> userDist(0, users.size() - 1);
        std::set<std::pair<int, int>> selectedPairs;
        std::set<int> usedUsers;

        while (selectedPairs.size() < 6) {
            int user1 = userDist();
            int user2 = userDist();

            // Убедимся, что оба пользователя не использованы и не равны друг другу
            if (user1 != user2 && usedUsers.find(user1) == usedUsers.end() && usedUsers.find(user2) == usedUsers.end()) {
                selectedPairs.emplace(std::min(user1, user2), std::max(user1, user2));
                usedUsers.insert(user1);
                usedUsers.insert(user2);
            }
        }

        for (const auto& pair : selectedPairs) {
            const UserTerminal& transmitter = users[pair.first];
            const UserTerminal& receiver = users[pair.second];

            std::cout << "// UT " << transmitter.id << " - UT " << receiver.id << " // " << std::endl << std::endl;

            // Расчёт лос угов и соовтетствующих значений ДН полей приёмника и передатчика
            double losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA;
            transmitter.calculateLOSAngles(transmitter, receiver, losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA);
            //std::cout << losPhiAOD << ", " << losThetaZOD << ", " << losPhiAOA << ", " << losThetaZOA << "\n ";

            //_______________________________________STEP_2________________________________________//
            // Вычисление расстояния между передатчиком и приемником
            double distance_tx_rx = calculateDistance(transmitter, receiver);
            //std::cout << "Distance between transmitter and receiver is " << distance_tx_rx << std::endl;

            //Определяем вид линка
            bool los = generateLOSorNLOS(distance_tx_rx,1);

            //__STEP3__//
            double path_loss = calculatePathLoss(los, distance_tx_rx, nu);
            //std::cout << "PathLoss[dB]  = " << path_loss << std::endl << std::endl;

            //_________________________________________________STEP_4__________________________________________//
            //std::cout << "\nLSP for LOS for User " << transmitter.id << " and User " << receiver.id << " :\n\n";
            LargeScaleParameters lsp(los, nu); //3ГГц
            //lsp.showParameters();
            /*
            if (los)
            {
                SF_LOS_CDF << lsp.shadowFading << std::endl;
                K_LOS_CDF << lsp.riceanK << std::endl;
                DS_LOS_CDF << lsp.delaySpread << std::endl;
                ASA_LOS_CDF << lsp.azimuthSpreadArrival << std::endl;
                ASD_LOS_CDF << lsp.azimuthSpreadDeparture << std::endl;
                ZSA_LOS_CDF << lsp.zenithSpreadArrival << std::endl;
                ZSD_LOS_CDF << lsp.zenithSpreadDeparture << std::endl;


            }
            else
            {
                SF_NLOS_CDF << lsp.shadowFading << std::endl;
                DS_NLOS_CDF << lsp.delaySpread << std::endl;
                ASA_NLOS_CDF << lsp.azimuthSpreadArrival << std::endl;
                ASD_NLOS_CDF << lsp.azimuthSpreadDeparture << std::endl;
                ZSA_NLOS_CDF << lsp.zenithSpreadArrival << std::endl;
                ZSD_NLOS_CDF << lsp.zenithSpreadDeparture << std::endl;


            }
            */



            //______________STEP_5_______________//
        // Генерация задержек кластеров
            std::vector<double> clusterDelays = generateClusterDelays(los, lsp.delaySpread, lsp.riceanK); // Передаем delaySpread из LSP

            //_____________STEP_6_______________//
            // Генерация мощностей кластеров
            std::vector<double> clusterPowers = generateClusterPowers(los, clusterDelays, lsp.delaySpread);

            std::vector<double> clusterPowersWithScalingFactors = clusterPowers;

            if (los) {

                for (size_t n = 0; n < clusterPowersWithScalingFactors.size(); ++n) {
                    clusterPowersWithScalingFactors[n] *= (1 / (pow(10, lsp.riceanK / 10) + 1));
                }
                clusterPowersWithScalingFactors[0] += pow(10, lsp.riceanK / 10) / (pow(10, lsp.riceanK / 10) + 1);



                std::vector<double> clusterPowers_main;

                double threshold = *max_element(clusterPowersWithScalingFactors.begin(), clusterPowersWithScalingFactors.end()) * 0.00316; // Порог для удаления
                for (size_t n = 0; n < clusterPowersWithScalingFactors.size(); ++n) {
                    if (clusterPowersWithScalingFactors[n] > threshold)
                    {
                        clusterPowers_main.emplace_back(clusterPowersWithScalingFactors[n]);

                    }
                    else { indicesToDelete.push_back(n); }
                }

                clusterPowersWithScalingFactors = clusterPowers_main;
            }


            //оставляю лишь нужные задержки  и сортирую их с мощностями 
            for (int i = indicesToDelete.size() - 1; i >= 0; --i) {
                clusterDelays.erase(clusterDelays.begin() + indicesToDelete[i]);
                clusterPowers.erase(clusterPowers.begin() + indicesToDelete[i]);
            }

            indicesToDelete.clear();

            

            //_____________STEP_7_______________//
            //Generate arrival angles and departure angles for both azimuth and elevation
            Eigen::MatrixXd PhiAOA = generateAOAorAOD_n_m(los, clusterPowersWithScalingFactors, lsp.azimuthSpreadArrival, lsp.riceanK, losPhiAOA, 0);
            Eigen::MatrixXd PhiAOD = generateAOAorAOD_n_m(los, clusterPowersWithScalingFactors, lsp.azimuthSpreadDeparture, lsp.riceanK, losPhiAOD, 1);
            Eigen::MatrixXd ThetaZOA = generateZOAorZOD_n_m(los, clusterPowersWithScalingFactors, lsp.zenithSpreadArrival, lsp.riceanK, losThetaZOA, 2);
            Eigen::MatrixXd ThetaZOD = generateZOAorZOD_n_m(los, clusterPowersWithScalingFactors, lsp.zenithSpreadDeparture, lsp.riceanK, losThetaZOD, 3);

            

            ///std::cout << "Coupling of rays within a cluster for both azimuth and elevation " << std::endl;
            randomCouplingRays(PhiAOD, PhiAOA, ThetaZOD, ThetaZOA, los);



            //___________A1___A2___________//
            std::vector<double> AS = calculateAngularSpreadandMeanAngles(los, clusterPowersWithScalingFactors, PhiAOD, PhiAOA, ThetaZOD, ThetaZOA);
            std::cout << "AS for AOD | AOA | ZOD | ZOA : ";
            for (size_t i = 0; i < 4 && i < AS.size(); ++i) {
                std::cout << AS[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Mean Angles for AOD | AOA | ZOD | ZOA : " << std::endl;
            for (size_t i = 4; i < 8 && i < AS.size(); ++i) {
                std::cout << AS[i] << " ";
            }
            std::cout << std::endl << std::endl;

            AS_AOA_CDF << AS[1] << std::endl;
            AS_AOD_CDF << AS[0] << std::endl;
            AS_ZOA_CDF << AS[3] << std::endl;
            AS_ZOD_CDF << AS[2] << std::endl;

            Mean_AOA_CDF << AS[5] << std::endl;
            Mean_AOD_CDF << AS[4] << std::endl;
            Mean_ZOA_CDF << AS[7] << std::endl;
            Mean_ZOD_CDF << AS[6] << std::endl;


        }
    }
    AS_AOA_CDF.close();
    AS_AOD_CDF.close();
    AS_ZOA_CDF.close();
    AS_ZOD_CDF.close();

    Mean_AOA_CDF.close();
    Mean_AOD_CDF.close();
    Mean_ZOA_CDF.close();
    Mean_ZOD_CDF.close();

    SF_LOS_CDF.close();
    K_LOS_CDF.close();
    DS_LOS_CDF.close();
    ASA_LOS_CDF.close();
    ASD_LOS_CDF.close();
    ZSA_LOS_CDF.close();
    ZSD_LOS_CDF.close();



    SF_NLOS_CDF.close();
    DS_NLOS_CDF.close();
    ASA_NLOS_CDF.close();
    ASD_NLOS_CDF.close();
    ZSA_NLOS_CDF.close();
    ZSD_NLOS_CDF.close();
    return 0;
}



/* построение на матлабе
filename = 'C:\Users\RadioChelik322\source\repos\Indoor_Office\SF_NLOS_CDF.txt';
data = load(filename);
sortedData = sort(data);
n = length(sortedData);
cdfValues = (1:n) / n;

figure;
plot(sortedData, cdfValues, 'LineWidth', 2);
xlabel('Значения');
ylabel('CDF');
title('Кумулятивная распределительная функция (CDF)');
grid on;
*/

/*
data1 = load('C:\Users\RadioChelik322\source\repos\Indoor_Office\DS_LOS_CDF.txt');
data2 = load('C:\Users\RadioChelik322\source\repos\Indoor_Office\K_LOS_CDF.txt');
correlationMatrix = corrcoef(data1, data2);
correlationCoefficient = correlationMatrix(1, 2);
disp(['Коэффициент корреляции: ', num2str(correlationCoefficient)]);
*/

/*
% Загрузка данных для LOS
filename_LOS = 'C:\Users\RadioChelik322\source\repos\Test_gen\AS_AOD_CDF.txt';
data_LOS = load(filename_LOS);
sortedData_LOS = sort(data_LOS);
n_LOS = length(sortedData_LOS);
cdfValues_LOS = (1:n_LOS) / n_LOS;

% Загрузка данных для NLOS
filename_NLOS = 'C:\Users\RadioChelik322\source\repos\Test_gen\AS_ZOD_CDF.txt';
data_NLOS = load(filename_NLOS);
sortedData_NLOS = sort(data_NLOS);
n_NLOS = length(sortedData_NLOS);
cdfValues_NLOS = (1:n_NLOS) / n_NLOS;

% Построение графиков
figure;
plot(sortedData_LOS, cdfValues_LOS, 'LineWidth', 2, 'DisplayName', 'AS_AOD'); % График LOS
hold on; % Удерживаем текущий график
plot(sortedData_NLOS, cdfValues_NLOS, 'LineWidth', 2, 'DisplayName', 'AS_ZOD'); % График NLOS

% Настройка графика
xlabel('Значения');
ylabel('CDF');
title('Кумулятивная распределительная функция (CDF)');
grid on;
legend show; % Показываем легенду
legend('AS_AOD', 'AS_ZOD', 'Location', 'Best');*/