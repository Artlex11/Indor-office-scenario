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

using namespace Eigen;
// Создание пары векторов

static std::pair<std::vector<double>, std::vector<double>> sort_with_indices(const std::vector<double>& vector1, const std::vector<double>& vector2)
{
    // Проверяем, что оба вектора имеют одинаковую длину
    if (vector1.size() != vector2.size()) {
        /*throw - обработка исключения */
        throw std::invalid_argument("Оба вектора должны иметь одинаковую длину.");
    }

    // Создаем вектор индексов
    std::vector<int> sorted_indices(vector1.size());
    for (int i = 0; i < vector1.size(); ++i) {
        sorted_indices[i] = i;
    }

    // Сортируем индексы по значениям первого вектора
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&vector1](int i1, int i2) {
        return vector1[i1] > vector1[i2]; // Сортировка по убыванию
        });

    // Создаем отсортированные векторы
    std::vector<double> sorted_vector1;
    std::vector<double> sorted_vector2;
    for (int index : sorted_indices) {
        sorted_vector1.push_back(vector1[index]);
        sorted_vector2.push_back(vector2[index]);
    }

    return { sorted_vector1, sorted_vector2 };
}

std::vector<int> indicesToDelete;

//______Функция_для_генерации_случайных_значений_x_по_нормальному_распределению___________//
double generateNormalRandom(double mean, double stddev, std::mt19937& gen)
{
    std::normal_distribution<double> dis(mean, stddev);
    return dis(gen); // Генерируем случайное значение по нормальному распределению
}

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
        std::random_device rd;
        std::mt19937 gen(rd());



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
            VectorXd standardDeviation(7);
            MatrixXd C(7, 7);

            //Если LOS :
            //Затемнение SF
            StandardDeviationSF = 3;
            MeanSF = 0;
            //std::normal_distribution<> loSshadowFadingDist(MeanSF, StandardDeviationSF);
            //shadowFading = loSshadowFadingDist(gen);

            //К-фактор (K)
            StandardDeviationK = 4;
            MeanK = 7;
            //std::normal_distribution<> riceanKDist0(MeanK, StandardDeviationK);
            //riceanK = (riceanKDist0(gen));

            // разброс задержки (DS)
            StandardDeviationDS = 0.18;
            MeanDS = (-0.01 * log10(1 + fc) - 7.692);
            //std::normal_distribution<> loSdelaySpreadDist(MeanDS, StandardDeviationDS);
            //delaySpread = loSdelaySpreadDist(gen);

            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.18;
            MeanASD = 1.6;
            //std::normal_distribution<> loSAzimuthSpreadDepartureDist(MeanASD, StandardDeviationASD);
            //azimuthSpreadDeparture = std::min((loSAzimuthSpreadDepartureDist(gen)), log10(104.0));

            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.119);
            MeanASA = (-0.19 * log10(1 + fc) + 1.781);
            //std::normal_distribution<> loSAzimuthSpreadArrivalDist(MeanASA, StandardDeviationASA);
            //azimuthSpreadArrival = std::min((loSAzimuthSpreadDepartureDist(gen)), log10(104.0));


            //For frequencies below 6 GHz, use fc = 6 when determining the values of the frequency-dependent
            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.04 * log10(1 + fc) + 0.264);
            MeanZSA = (-0.26 * log10(1 + fc) + 1.44);
            //std::normal_distribution<> loSZenithSpreadArrivalDist(MeanZSA, StandardDeviationZSA);
            //zenithSpreadArrival = std::min(loSZenithSpreadArrivalDist(gen), log10(52.0));

            //Зенитный угол Разброс вылета (ZSD)

            if (fc < 6.0) {
                StandardDeviationZSD = (0.13 * log10(1 + 6) + 0.30);
                MeanZSD = (-1.43 * log10(1 + 6) + 2.228);
                //std::normal_distribution<> loSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
                //zenithSpreadDeparture = std::min(loSZenithSpreadDepartureDist(gen), log10(52.0));
            }
            else {
                StandardDeviationZSD = (0.13 * log10(1 + fc) + 0.30);
                MeanZSD = (-1.43 * log10(1 + fc) + 2.228);
                //std::normal_distribution<> loSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
                //zenithSpreadDeparture = std::min(loSZenithSpreadDepartureDist(gen), log10(52.0));
            }
            /*
            shadowFading = pow(10, shadowFading / 10);
            riceanK = pow(10, riceanK / 10);
            delaySpread = pow(10, delaySpread);
            azimuthSpreadDeparture = pow(10, azimuthSpreadDeparture);
            azimuthSpreadArrival = pow(10, azimuthSpreadArrival);
            zenithSpreadDeparture = pow(10, zenithSpreadDeparture);
            zenithSpreadArrival = pow(10, zenithSpreadArrival);
            */

            value << shadowFading, riceanK, delaySpread, azimuthSpreadDeparture, azimuthSpreadArrival, zenithSpreadDeparture, zenithSpreadArrival;
            means << MeanSF, MeanK, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;
            standardDeviation << StandardDeviationSF, StandardDeviationK, StandardDeviationDS, StandardDeviationASD, StandardDeviationASA, StandardDeviationZSD, StandardDeviationZSA;
            
            

            //SF K DS ASD ASA ZSD ZSA
            
            C << StandardDeviationSF * StandardDeviationSF, 0.5 * StandardDeviationSF * StandardDeviationK, -0.8 * StandardDeviationSF * StandardDeviationDS, -0.4 * StandardDeviationSF * StandardDeviationASD, -0.5 * StandardDeviationSF * StandardDeviationASA, 0.2 * StandardDeviationSF * StandardDeviationZSD, 0.3 * StandardDeviationSF * StandardDeviationZSA,
                0.5 * StandardDeviationK * StandardDeviationSF, StandardDeviationK* StandardDeviationK, -0.5 * StandardDeviationK * StandardDeviationDS, 0.0 * StandardDeviationK * StandardDeviationASD, 0.0 * StandardDeviationK * StandardDeviationASA, 0.0 * StandardDeviationK * StandardDeviationZSD, 0.1 * StandardDeviationK * StandardDeviationZSA,
                -0.8 * StandardDeviationDS * StandardDeviationSF, -0.5 * StandardDeviationDS * StandardDeviationK, StandardDeviationDS* StandardDeviationDS, 0.6 * StandardDeviationDS * StandardDeviationASD, 0.8 * StandardDeviationDS * StandardDeviationASA, 0.1 * StandardDeviationDS * StandardDeviationZSD, 0.2 * StandardDeviationDS * StandardDeviationZSA,
                -0.4 * StandardDeviationASD * StandardDeviationSF, 0.0 * StandardDeviationASD * StandardDeviationK, 0.6 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.4 * StandardDeviationASD * StandardDeviationASA, 0.5 * StandardDeviationASD * StandardDeviationZSD, 0.0 * StandardDeviationASD * StandardDeviationZSA,
                -0.5 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationK, 0.8 * StandardDeviationASA * StandardDeviationDS, 0.4 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, 0.0 * StandardDeviationASA * StandardDeviationZSD, 0.5 * StandardDeviationASA * StandardDeviationZSA,
                0.2 * StandardDeviationZSD * StandardDeviationSF, 0.0 * StandardDeviationZSD * StandardDeviationK, 0.1 * StandardDeviationZSD * StandardDeviationDS, 0.5 * StandardDeviationZSD * StandardDeviationASD, 0.0 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.0 * StandardDeviationZSD * StandardDeviationZSA,
                0.3 * StandardDeviationZSA * StandardDeviationSF, 0.1 * StandardDeviationZSA * StandardDeviationK, 0.2 * StandardDeviationZSA * StandardDeviationDS, 0.0 * StandardDeviationZSA * StandardDeviationASD, 0.5 * StandardDeviationZSA * StandardDeviationASA, 0.0 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;
            /*
            C << 1, 0.5 , -0.8 , -0.4 , -0.5, 0.2 , 0.3 ,
                0.5 , 1, -0.5 , 0.0 , 0.0 , 0.0, 0.1 ,
                -0.8 , -0.5 , 1, 0.6 , 0.8 , 0.1 , 0.2 ,
                -0.4 , 0.0, 0.6 , 1, 0.4 , 0.5 , 0.0 ,
                -0.5 , 0.0 , 0.8 , 0.4 , 1, 0.0 , 0.5 ,
                0.2 , 0.0 , 0.1 , 0.5 , 0.0 , 1, 0.0 ,
                0.3 , 0.1 , 0.2 , 0.0 , 0.5 , 0.0 , 1;
            */

            //std::cout << C << std::endl << std::endl;


            /*Вычисление и проверка главных миноров
            for (int i = 1; i <= C.rows(); ++i) {
                // Получаем верхнюю левую подматрицу размером i x i
                MatrixXd minor = C.topLeftCorner(i, i);

                // Вычисляем определитель (главный минор)
                double determinant = minor.determinant();

                // Выводим значение и знак главного минора
                std::cout << i << "x" << i << ": " << determinant << std::endl;
                if (determinant > 0) {
                    std::cout << "+ " << std::endl;
                }
                else if (determinant < 0) {
                    std::cout << "- " << std::endl;
                }
                else {
                    std::cout << "0 " << std::endl;
                }
            }
            */

            MatrixXd L;
            L = C.llt().matrixL();
            //std::cout << L << "\n";

            // Проверка на положительную определенность
            if (C.llt().info() != Success) {
                std::cerr << "Matrix C is not positive definite!" << std::endl;
                return;
            }

            
            std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
            std::normal_distribution<double> normDist(0.0, 1.0);
            for (int i = 0; i < 7; ++i) {
                value(i) = normDist(gen);
            }

            // Преобразование случайных величин с помощью матрицы Лапласа

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
            VectorXd standardDeviation(6);
            MatrixXd C(6, 6);

            //Если NLOS:
            //Затенение (SF)
            StandardDeviationSF = 8.03;
            MeanSF = 0;
            //std::normal_distribution<> nLoSshadowFadingDist(MeanSF, StandardDeviationSF);
            //shadowFading = nLoSshadowFadingDist(gen);

            // разброс задержки (DS)
            StandardDeviationDS = (0.10 * log10(1 + fc) + 0.055);
            MeanDS = (-0.28 * log10(1 + fc) - 7.173);
            //std::normal_distribution<> nLoSdelaySpreadDist(MeanDS, StandardDeviationDS);
            //delaySpread = nLoSdelaySpreadDist(gen); // Задержка распространения НЛОС

            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.25;
            MeanASD = 1.62;
            //std::normal_distribution<> nLoSAzimuthSpreadDepartureDist(MeanASD, StandardDeviationASD);
            //azimuthSpreadDeparture = std::min((nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));

            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.059);
            MeanASA = (-0.11 * log10(1 + fc) + 1.863);
            //std::normal_distribution<> nLoSAzimuthSpreadArrivaDist(MeanASA, StandardDeviationASA);
            //azimuthSpreadArrival = std::min((nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));


            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.09 * log10(1 + fc) + 0.746);
            MeanZSA = (-0.15 * log10(1 + fc) + 1.387);
            //std::normal_distribution<> nLoSZenithSpreadArrivalDist(MeanZSA, StandardDeviationZSA);
            //zenithSpreadArrival = std::min(nLoSZenithSpreadArrivalDist(gen), log10(52.0));

            //Зенитный угол Разброс вылета (ZSD)
            StandardDeviationZSD = 0.36;
            MeanZSD = 1.08;

            //std::normal_distribution<> nLoSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
            //zenithSpreadDeparture = std::min(nLoSZenithSpreadDepartureDist(gen), log10(52.0));
            
            value << shadowFading, delaySpread, azimuthSpreadDeparture, azimuthSpreadArrival, zenithSpreadDeparture, zenithSpreadArrival;
            means << MeanSF, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;
            standardDeviation << StandardDeviationSF,  StandardDeviationDS, StandardDeviationASD, StandardDeviationASA, StandardDeviationZSD, StandardDeviationZSA;

            std::cout << value << "\n";
            

            //SF  DS ASD ASA ZSD ZSA
            C << StandardDeviationSF * StandardDeviationSF, -0.5 * StandardDeviationSF * StandardDeviationDS, 0.0 * StandardDeviationSF * StandardDeviationASD, -0.4 * StandardDeviationSF * StandardDeviationASA, 0.0 * StandardDeviationSF * StandardDeviationZSD, 0.0 * StandardDeviationSF * StandardDeviationZSA,
                -0.5 * StandardDeviationDS * StandardDeviationSF, StandardDeviationDS* StandardDeviationDS, 0.4 * StandardDeviationDS * StandardDeviationASD, 0.0 * StandardDeviationDS * StandardDeviationASA, -0.27 * StandardDeviationDS * StandardDeviationZSD, -0.06 * StandardDeviationDS * StandardDeviationZSA,
                0.0 * StandardDeviationASD * StandardDeviationSF, 0.4 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.0 * StandardDeviationASD * StandardDeviationASA, 0.35 * StandardDeviationASD * StandardDeviationZSD, 0.23 * StandardDeviationASD * StandardDeviationZSA,
                -0.4 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationDS, 0.0 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, -0.08 * StandardDeviationASA * StandardDeviationZSD, 0.43 * StandardDeviationASA * StandardDeviationZSA,
                0.0 * StandardDeviationZSD * StandardDeviationSF, -0.27 * StandardDeviationZSD * StandardDeviationDS, 0.35 * StandardDeviationZSD * StandardDeviationASD, -0.08 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.42 * StandardDeviationZSD * StandardDeviationZSA,
                0.0 * StandardDeviationZSA * StandardDeviationSF, -0.06 * StandardDeviationZSA * StandardDeviationDS, 0.23 * StandardDeviationZSA * StandardDeviationASD, 0.43 * StandardDeviationZSA * StandardDeviationASA, 0.42 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;
            
            /*
            C << 1, -0.5 , 0.0 , -0.4 , 0.0 , 0.0 ,
                -0.5 , 1, 0.4 , 0.0 , -0.27 , -0.06 ,
                0.0 , 0.4 , 1, 0.0 , 0.35 , 0.23 ,
                -0.4 , 0.0 , 0.0 , 1, -0.08 , 0.43 ,
                0.0 , -0.27 , 0.35 , -0.08 , 1 , 0.42 ,
                0.0 , -0.06 , 0.23 , 0.43 , 0.42 , 1;
            */
            //std::cout << C << std::endl << std::endl;


            /* Вычисление и проверка главных миноров
            for (int i = 1; i <= C.rows(); ++i) {
                // Получаем верхнюю левую подматрицу размером i x i
                MatrixXd minor = C.topLeftCorner(i, i);

                // Вычисляем определитель (главный минор)
                double determinant = minor.determinant();

                // Выводим значение и знак главного минора
                std::cout << i << "x" << i << ": " << determinant << std::endl;
                if (determinant > 0) {
                    std::cout << "+ " << std::endl;
                }
                else if (determinant < 0) {
                    std::cout << "- " << std::endl;
                }
                else {
                    std::cout << "0 " << std::endl;
                }
            }
            */

            MatrixXd L;
            L = C.llt().matrixL();
            //std::cout << L << "\n";

            // Проверка на положительную определенность
            if (C.llt().info() != Success) {
                std::cerr << "Matrix C is not positive definite!" << std::endl;
                return;
            }

            
            std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
            std::normal_distribution<double> normDist(0.0, 1.0);
            for (int i = 0; i < 6; ++i) {
                value(i) = normDist(gen);
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



        double A3D_PowerPattern = 8 - std::min((-1) * (AverticalPowerPattern + AhorizontalPowerPattern), 30.0); // знак "-" в начале формулы

        A3D_PowerPattern = pow(10, A3D_PowerPattern / 20);


        double f1_Theta = sqrt(A3D_PowerPattern) * cos(ksi);
        double f1_Phi = sqrt(A3D_PowerPattern) * sin(ksi);

        fieldPattern << f1_Theta, f1_Phi;

        return fieldPattern;
    }

    // в fieldPattern - 1 theta, 2 phi
    /*
    //Переход от LSC в GSC
    Eigen::Vector2d transformationFromLCStoGCS(double thetaAngle, double phiAngle, double bearingAngle, double downtiltAngle, double slantAngle, Eigen::Vector2d& fieldPattern) const
    {
        Eigen::Vector2d transformFieldPattern;

        double cos_Pci = (cos(downtiltAngle) * cos(slantAngle) * sin(thetaAngle) - (sin(downtiltAngle) * cos(slantAngle) * cos(phiAngle - bearingAngle) - sin(slantAngle) * sin(phiAngle - bearingAngle)) * cos(thetaAngle))
            / sqrt(1 - pow(cos(downtiltAngle) * cos(slantAngle) * sin(thetaAngle) + (sin(downtiltAngle) * cos(slantAngle) * cos(phiAngle - bearingAngle) - sin(slantAngle) * sin(phiAngle - bearingAngle)) * sin(thetaAngle), 2));


        double sin_Pci = (sin(downtiltAngle) * cos(slantAngle) * sin(phiAngle - bearingAngle) + sin(slantAngle) * cos(phiAngle - bearingAngle))
            / sqrt(1 - pow(cos(downtiltAngle) * cos(slantAngle) * sin(thetaAngle) + (sin(downtiltAngle) * cos(slantAngle) * cos(phiAngle - bearingAngle) - sin(slantAngle) * sin(phiAngle - bearingAngle)) * sin(thetaAngle),2));

        double f_Theta = cos_Pci * fieldPattern(0) - sin_Pci * fieldPattern(1);
        double f_Phi = sin_Pci * fieldPattern(0) + cos_Pci * fieldPattern(1);
        //double f_Theta = sin_Pci * fieldPattern(0) + cos_Pci * fieldPattern(1);
        //double f_Phi = cos_Pci * fieldPattern(0) - sin_Pci * fieldPattern(1);

        transformFieldPattern << f_Theta, f_Phi;
        return transformFieldPattern;
    }
    */
    void calculateLOSAngles(const UserTerminal& transmitter, const UserTerminal& receiver,
        double& losPhiAOD, double& losThetaZOD,
        double& losPhiAOA, double& losThetaZOA) const
    {
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

    Eigen::MatrixXd generateAntennaElements() const
    {
        Eigen::MatrixXd locationMatrixAntennaElements(16, 3);

        locationMatrixAntennaElements.row(0) << -3 * wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, -3 * wavelength * cos(bearingAngle) * cos(slantAngle) / 4, wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(8) = locationMatrixAntennaElements.row(0);
        locationMatrixAntennaElements.row(1) << -wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, -wavelength * cos(bearingAngle) * cos(slantAngle) / 4, wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(9) = locationMatrixAntennaElements.row(1);
        locationMatrixAntennaElements.row(2) << wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, wavelength* cos(bearingAngle)* cos(slantAngle) / 4, wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(10) = locationMatrixAntennaElements.row(2);
        locationMatrixAntennaElements.row(3) << 3 * wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, 3 * wavelength * cos(bearingAngle) * cos(slantAngle) / 4, wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(11) = locationMatrixAntennaElements.row(3);
        locationMatrixAntennaElements.row(4) << -3 * wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, -3 * wavelength * cos(bearingAngle) * cos(slantAngle) / 4, -wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(12) = locationMatrixAntennaElements.row(4);
        locationMatrixAntennaElements.row(5) << -wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, -wavelength * cos(bearingAngle) * cos(slantAngle) / 4, -wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(13) = locationMatrixAntennaElements.row(5);
        locationMatrixAntennaElements.row(6) << wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, wavelength* cos(bearingAngle)* cos(slantAngle) / 4, -wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(14) = locationMatrixAntennaElements.row(6);
        locationMatrixAntennaElements.row(7) << 3 * wavelength * cos(M_PI / 2 + bearingAngle) * cos(downtiltAngle) / 4, 3 * wavelength * cos(bearingAngle) * cos(slantAngle) / 4, -wavelength / 4 * cos(downtiltAngle) * cos(slantAngle);
        locationMatrixAntennaElements.row(15) = locationMatrixAntennaElements.row(7);
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

std::vector<double> generateClusterDelays(bool los, double delaySpread,  double riceanK)
{
    delaySpread = pow(10, delaySpread);

    

    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
    

    std::vector<double> delays_tau;
    double delay_tau_n;
    //delays_tau.clear();

    if (los) {
        for (int n = 0; n < 15; ++n)
        {
            double Xn = std::uniform_real_distribution<>(0.0, 1.0)(gen); // Xn ~ uniform(0,1)
            delay_tau_n = -1 * los_r_tau * log(Xn) * delaySpread;
            delays_tau.push_back(delay_tau_n);
        }
    }
    else {
        for (int n = 0; n < 19; ++n)
        {
            double Xn = std::uniform_real_distribution<>(0.0, 1.0)(gen); // Xn ~ uniform(0,1)
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
    /*
    // Если LOS, применяем масштабирование (The scaled delays are not to be used in cluster power generation. )
    if (riceanK > 0) { // Если K-фактор положителен
        double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) + (0.0433 *riceanK) + 0.7705;

        delay[0] /= scalingFactor; // Масштабирование задержек
    }
    */
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
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());

    for (size_t n = 0; n < clusterDelays.size(); ++n) {

        if (los) {
            std::normal_distribution<> shadowFadingDist(0, 36); // is the per cluster shadowing term in [dB] LOS.
            double shadowing = shadowFadingDist(gen); // Генерация затенения
            power = exp(((-1) * clusterDelays[n] * (los_r_tau - 1))  * pow(10, ( (-1) * shadowing / 10)) / (los_r_tau * delaySpread));
        }
        else {
            std::normal_distribution<> shadowFadingDist(0, 9); // is the per cluster shadowing term in [dB] NLOS.
            double shadowing = shadowFadingDist(gen); // Генерация затенения
            power = exp(((-1) * clusterDelays[n] * (nlos_r_tau - 1))  * pow(10, ( (-1) * shadowing / 10)) / (nlos_r_tau * delaySpread));
        }
        clusterPowers[n] = power;
    }


    // Нормализация мощностей кластеров
    
    for (auto& n : clusterPowers) sumclusterPowers += n;

    
    for (size_t n = 0; n < clusterPowers.size(); ++n) {
        clusterPowers[n] = clusterPowers[n] / sumclusterPowers; // Нормализуем по суммарной мощности    
    }
    maxPower = *max_element(clusterPowers.begin(), clusterPowers.end());

    std::vector<double> clusterPowers_main;

    double threshold = maxPower * 0.00316; // Порог для удаления
    for (size_t n = 0; n < clusterPowers.size(); ++n) {
        if (clusterPowers[n] > threshold)
        {
            clusterPowers_main.emplace_back(clusterPowers[n]);

        }
        else { indicesToDelete.push_back(n); }
    }

    return clusterPowers_main; // Возвращаем нормализованные мощности кластеров
};

//________________________________________________STEP_7___________________________________________________//
//____________________________Cluster_ASA_,_ASD_,_ZSA_in_[deg]_for_los_and_nlos____________________________//
double los_C_ASA = 8.0;
double los_C_ASD = 5.0;
double los_C_ZSA = 9.0;
double los_C_ZSD = 0.375;

double nlos_C_ASA = 11.0;
double nlos_C_ASD = 5.0;
double nlos_C_ZSA = 9.0;
double nlos_C_ZSD = 0.375;

//____table_7.5.3________Ray_offset_angles_within_a_cluster,_given_for_rms_angle_spread_normalized_to_1
std::vector<double> am = { 0.0447, -0.0447 ,0.1413 ,-0.1413,0.2492 ,-0.2492 ,0.3715 ,-0.3715 ,0.5129 ,-0.5129 ,
    0.6797 ,-0.6797, 0.8844 , -0.8844,1.1481,-1.1481, 1.5195 , -1.5195, 2.1551 , -2.1551 };

//_________________________________________________AOA____________________________________________________//
Eigen::MatrixXd generatePhiAOA(bool los, const std::vector<double>& clusterPowers_main, double AzimuthSpreadArrival, double riceanK, double losPhiAOA) {

    losPhiAOA = losPhiAOA * 180 / M_PI;
    AzimuthSpreadArrival = pow(10, AzimuthSpreadArrival);
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());

    double maxPower = *max_element(clusterPowers_main.begin(), clusterPowers_main.end());
    //n - должен быть размер кластера задержек
    double X1;
    double Y1;
    double Phi_1_AOA;

    double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

    if (los) // случай LOS
    {
        std::vector<double> Phi_n_AOA(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_AOA(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            if (los)
            {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                double Yn = generateNormalRandom(0, (AzimuthSpreadArrival / 7) * (AzimuthSpreadArrival / 7), gen);
                // Yn ~ N(0,(ASA/7)^2)

                if (n == 0) // добавляем 1 луч
                {
                    X1 = Xn; Y1 = Yn; Phi_1_AOA = Phi_n_AOA[n];
                    Phi_n_m_AOA(n, 0) = losPhiAOA;
                }

                else
                {
                    Phi_n_AOA[n] = (2 * (AzimuthSpreadArrival / 1.4) * pow(-log(clusterPowers_main[n - 1] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 

                    Phi_n_AOA[n] = Phi_n_AOA[n] * Xn + Yn + losPhiAOA - Phi_1_AOA * X1 - Y1;
                    for (int m = 0; m < 20; ++m)
                    {
                        Phi_n_m_AOA(n, m) = Phi_n_AOA[n] + los_C_ASA * am[m];
                    }
                }
            }
        }
        return Phi_n_m_AOA;
    }
    else // случай NLOS
    {
        std::vector<double> Phi_n_AOA(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_AOA(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            double C_phi = 1.273;
            Phi_n_AOA[n] = (2 * (AzimuthSpreadArrival / 1.4) * pow(-log(clusterPowers_main[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
            double Yn = generateNormalRandom(0, (AzimuthSpreadArrival / 7) * (AzimuthSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

            Phi_n_AOA[n] = Phi_n_AOA[n] * Xn + Yn + losPhiAOA;

            for (int m = 0; m < 20; ++m)
            {
                Phi_n_m_AOA(n, m) = Phi_n_AOA[n] + nlos_C_ASA * am[m];
            }
        }
        return Phi_n_m_AOA;
    }
}

//_____________________________________________AOD____________________________________________________//
Eigen::MatrixXd generatePhiAOD(bool los, const std::vector<double>& clusterPowers_main, double AzimuthSpreadDeparture, double riceanK, double losPhiAOD) {

    losPhiAOD = losPhiAOD * 180 / M_PI;
    AzimuthSpreadDeparture = pow(10, AzimuthSpreadDeparture);
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());

    double maxPower = *max_element(clusterPowers_main.begin(), clusterPowers_main.end());
    //n - должен быть размер кластера задержек
    double X1;
    double Y1;
    double Phi_1_AOD;

    double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

    if (los) // случай LOS
    {
        std::vector<double> Phi_n_AOD(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_AOD(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            if (los)
            {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                double Yn = generateNormalRandom(0, (AzimuthSpreadDeparture / 7) * (AzimuthSpreadDeparture / 7), gen);
                // Yn ~ N(0,(ASD/7)^2)

                if (n == 0) // добавляем 1 луч
                {
                    X1 = Xn; Y1 = Yn; Phi_1_AOD = Phi_n_AOD[n];
                    Phi_n_m_AOD(n, 0) = losPhiAOD;
                }

                else
                {
                    Phi_n_AOD[n] = (2 * (AzimuthSpreadDeparture / 1.4) * pow(-log(clusterPowers_main[n - 1] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 

                    Phi_n_AOD[n] = Phi_n_AOD[n] * Xn + Yn + losPhiAOD - Phi_1_AOD * X1 - Y1;
                    for (int m = 0; m < 20; ++m)
                    {
                        Phi_n_m_AOD(n, m) = Phi_n_AOD[n] + los_C_ASD * am[m];
                    }
                }
            }
        }
        return Phi_n_m_AOD;
    }
    else // случай NLOS
    {
        std::vector<double> Phi_n_AOD(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_AOD(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            double C_phi = 1.273;
            Phi_n_AOD[n] = (2 * (AzimuthSpreadDeparture / 1.4) * pow(-log(clusterPowers_main[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
            double Yn = generateNormalRandom(0, (AzimuthSpreadDeparture / 7) * (AzimuthSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

            Phi_n_AOD[n] = Phi_n_AOD[n] * Xn + Yn + losPhiAOD;

            for (int m = 0; m < 20; ++m)
            {
                Phi_n_m_AOD(n, m) = Phi_n_AOD[n] + nlos_C_ASD * am[m];
            }
        }
        return Phi_n_m_AOD;
    }
}

//_____________________________________________ZOA____________________________________________________//
Eigen::MatrixXd generateThetaZOA(bool los, const std::vector<double>& clusterPowers_main, double ZenithSpreadArrival, double riceanK, double losPhiZOA) {

    losPhiZOA = losPhiZOA * 180 / M_PI;
    ZenithSpreadArrival = pow(10, ZenithSpreadArrival);
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());

    double maxPower = *max_element(clusterPowers_main.begin(), clusterPowers_main.end());
    //n - должен быть размер кластера задержек
    double X1;
    double Y1;
    double Phi_1_ZOA;

    double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

    if (los) // случай LOS
    {
        std::vector<double> Phi_n_ZOA(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_ZOA(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            if (los)
            {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                double Yn = generateNormalRandom(0, (ZenithSpreadArrival / 7) * (ZenithSpreadArrival / 7), gen);
                // Yn ~ N(0,(ZSA/7)^2)

                if (n == 0) // добавляем 1 луч
                {
                    X1 = Xn; Y1 = Yn; Phi_1_ZOA = Phi_n_ZOA[n];
                    Phi_n_m_ZOA(n, 0) = losPhiZOA;
                }

                else
                {
                    Phi_n_ZOA[n] = (2 * (ZenithSpreadArrival / 1.4) * pow(-log(clusterPowers_main[n - 1] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 

                    Phi_n_ZOA[n] = Phi_n_ZOA[n] * Xn + Yn + losPhiZOA - Phi_1_ZOA * X1 - Y1;
                    for (int m = 0; m < 20; ++m)
                    {
                        Phi_n_m_ZOA(n, m) = Phi_n_ZOA[n] + los_C_ZSA * am[m];
                    }
                }
            }
        }
        return Phi_n_m_ZOA;
    }
    else // случай NLOS
    {
        std::vector<double> Phi_n_ZOA(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_ZOA(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            double C_phi = 1.273;
            Phi_n_ZOA[n] = (2 * (ZenithSpreadArrival / 1.4) * pow(-log(clusterPowers_main[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
            double Yn = generateNormalRandom(0, (ZenithSpreadArrival / 7) * (ZenithSpreadArrival / 7), gen);//Yn ~ N(0,(ASA/7)^2)

            Phi_n_ZOA[n] = Phi_n_ZOA[n] * Xn + Yn + losPhiZOA;

            for (int m = 0; m < 20; ++m)
            {
                Phi_n_m_ZOA(n, m) = Phi_n_ZOA[n] + nlos_C_ZSA * am[m];
            }
        }
        return Phi_n_m_ZOA;
    }
}

//_________________________________________________ZOD____________________________________________________//
Eigen::MatrixXd generateThetaZOD(bool los, const std::vector<double>& clusterPowers_main, double ZenithSpreadDeparture, double riceanK, double losPhiZOD) {

    losPhiZOD = losPhiZOD * 180 / M_PI;
    ZenithSpreadDeparture = pow(10, ZenithSpreadDeparture);
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());

    double maxPower = *max_element(clusterPowers_main.begin(), clusterPowers_main.end());
    //n - должен быть размер кластера задержек
    double X1;
    double Y1;
    double Phi_1_ZOD;

    double Xn = std::uniform_real_distribution<>(-1.0, 1.0)(gen); // Xn ~ uniform(-1,1)

    if (los) // случай LOS
    {
        std::vector<double> Phi_n_ZOD(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_ZOD(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            if (los)
            {
                double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) + (0.028 * riceanK) + 1.035);
                double Yn = generateNormalRandom(0, (ZenithSpreadDeparture / 7) * (ZenithSpreadDeparture / 7), gen);
                // Yn ~ N(0,(ZSD/7)^2)

                if (n == 0) // добавляем 1 луч
                {
                    X1 = Xn; Y1 = Yn; Phi_1_ZOD = Phi_n_ZOD[n];
                    Phi_n_m_ZOD(n, 0) = losPhiZOD;
                }

                else
                {
                    Phi_n_ZOD[n] = (2 * (ZenithSpreadDeparture / 1.4) * pow(-log(clusterPowers_main[n - 1] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 

                    Phi_n_ZOD[n] = Phi_n_ZOD[n] * Xn + Yn + losPhiZOD - Phi_1_ZOD * X1 - Y1;
                    for (int m = 0; m < 20; ++m)
                    {
                        Phi_n_m_ZOD(n, m) = Phi_n_ZOD[n] + los_C_ZSD * am[m];
                    }
                }
            }
        }
        return Phi_n_m_ZOD;
    }
    else // случай NLOS
    {
        std::vector<double> Phi_n_ZOD(clusterPowers_main.size());
        Eigen::MatrixXd Phi_n_m_ZOD(clusterPowers_main.size(), 20);
        for (int n = 0; n < clusterPowers_main.size(); n++)
        {
            double C_phi = 1.273;
            Phi_n_ZOD[n] = (2 * (ZenithSpreadDeparture / 1.4) * pow(-log(clusterPowers_main[n] / maxPower), 0.5)) / C_phi; // Применение уравнения (7.5-9) 
            double Yn = generateNormalRandom(0, (ZenithSpreadDeparture / 7) * (ZenithSpreadDeparture / 7), gen);//Yn ~ N(0,(ASA/7)^2)

            Phi_n_ZOD[n] = Phi_n_ZOD[n] * Xn + Yn + losPhiZOD;

            for (int m = 0; m < 20; ++m)
            {
                Phi_n_m_ZOD(n, m) = Phi_n_ZOD[n] + nlos_C_ZSD * am[m];
            }
        }
        return Phi_n_m_ZOD;
    }
}
//_________________________________________A_1________________________________________________________//
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
        if (n == 0 && los) { continue; }
        for (int m = 0; m < 20; ++m) {
            weighted_sumAOD += (clasterPowers[n] / 20) * std::exp(std::complex<double>(0.0, AOD(n, m)));
            weighted_sumAOA += (clasterPowers[n] / 20) * std::exp(std::complex<double>(0.0, AOA(n, m)));
            weighted_sumZOD += (clasterPowers[n] / 20) * std::exp(std::complex<double>(0.0, ZOD(n, m)));
            weighted_sumZOA += (clasterPowers[n] / 20) * std::exp(std::complex<double>(0.0, ZOA(n, m)));
            weighted_sumPowers += (clasterPowers[n] / 20);

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

//______________________________________________STEP_9_____________________________________________________//
//___________________________________Функция_для_генерации_матриц_поляризации______________________________//
Eigen::MatrixXd generateXPR(bool los, const std::vector<double>& clusterPowers) {


    std::random_device rd;
    std::mt19937 gen(rd());

    Eigen::MatrixXd XPR(clusterPowers.size(), 20);
    double mean_Xn = 10;

    if (los) { mean_Xn = 11; }

    for (int n = 0; n < clusterPowers.size(); ++n) {
        for (int m = 0; m < 20; ++m) {

            std::normal_distribution<> X_n_m_Dist(mean_Xn, 16);
            XPR(n, m) = pow(10, (X_n_m_Dist(gen)) / 10);
        }
    }
    return XPR;
};

//______________________________________________STEP_10____________________________________________________//
//_________________________________Функция_для_генерации_случайных_начальных_фаз___________________________//
Eigen::MatrixXd generateInitialRandomPhases(std::vector<double>& clusterPowers)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> phaseDist(-M_PI, M_PI); // Распределение фаз от -π до π

    Eigen::MatrixXd initialRandomPhases(clusterPowers.size(), 80);

    for (int n = 0; n < clusterPowers.size(); ++n) {
        for (int m = 0; m < 20; ++m) {
            // Генерация случайных фаз для каждой комбинации поляризации
            initialRandomPhases(n, m * 4) = phaseDist(gen); // Theta - Theta
            initialRandomPhases(n, m * 4 + 1) = phaseDist(gen); // Theta - Phi
            initialRandomPhases(n, m * 4 + 2) = phaseDist(gen); // Phi - Theta
            initialRandomPhases(n, m * 4 + 3) = phaseDist(gen); // Phi - Phi
        }
    }
    return initialRandomPhases;
}

//______________________________________________STEP_11____________________________________________________//
//ИХ для НЛОС
Eigen::MatrixXcd generateNLOSChannelCoefficients(const UserTerminal& transmitter, const UserTerminal& receiver,
    std::vector<double>& clusterPowers, Eigen::MatrixXd& phiAOD_n_m, Eigen::MatrixXd& phiAOA_n_m,
    Eigen::MatrixXd& thetaZOD_n_m, Eigen::MatrixXd& thetaZOA_n_m, Eigen::MatrixXd& XRP, Eigen::MatrixXd& initialPhases) {

    phiAOD_n_m = phiAOD_n_m * M_PI / 180;
    phiAOA_n_m = phiAOA_n_m * M_PI / 180;
    thetaZOD_n_m = thetaZOD_n_m * M_PI / 180;
    thetaZOA_n_m = thetaZOA_n_m * M_PI / 180;

    std::complex<double> j(0.0, 1.0);
    Eigen::MatrixXcd channelCoefficients_u_s_n(256, clusterPowers.size() + 4);//16*16 комбинаций t-u
    channelCoefficients_u_s_n.setZero();

    // пары U-S
    for (int n = 0; n < clusterPowers.size(); ++n) {
        double ksi_rx;
        double ksi_tx;
        int pair = 0;
        for (int u = 0; u < 16; ++u) {
            if (u < 8) { ksi_rx = 45; }
            else { ksi_rx = -45; }
            for (int s = 0; s < 16; ++s) {
                if (s < 8) { ksi_tx = 45; }
                else { ksi_tx = -45; }
                for (int m = 0; m < 20; ++m) {
                    if (n > 1) {
                        Eigen::Vector2d F_tx = transmitter.FieldPattern(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), ksi_tx);
                        Eigen::Vector2d F_rx = receiver.FieldPattern(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), ksi_rx);

                        //Eigen::Vector2d F_tx = transmitter.transformationFromLCStoGCS(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F1_tx);
                        //Eigen::Vector2d F_rx = receiver.transformationFromLCStoGCS(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F1_rx);

                        Eigen::Vector3d sphericalUnitVector_tx(sin(thetaZOD_n_m(n, m)) * cos(phiAOD_n_m(n, m)),
                            sin(thetaZOD_n_m(n, m)) * sin(phiAOD_n_m(n, m)),
                            cos(thetaZOD_n_m(n, m)));
                        Eigen::Vector3d sphericalUnitVector_rx(sin(thetaZOA_n_m(n, m)) * cos(phiAOA_n_m(n, m)),
                            sin(thetaZOA_n_m(n, m)) * sin(phiAOA_n_m(n, m)),
                            cos(thetaZOA_n_m(n, m)));

                        Eigen::Matrix2cd XPR_and_InitialRandomPhases;
                        XPR_and_InitialRandomPhases <<
                            exp(j * initialPhases(n, m * 4)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 1)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 2)),
                            exp(j * initialPhases(n, m * 4 + 3));

                        double tx = sphericalUnitVector_tx.transpose() * transmitter.generateAntennaElements().row(s).transpose();
                        double rx = sphericalUnitVector_rx.transpose() * receiver.generateAntennaElements().row(u).transpose();

                        // Разделение сложной операции на более простые
                        auto temp1 = F_rx.transpose() * XPR_and_InitialRandomPhases; // Промежуточный результат
                        auto temp2 = temp1 * F_tx; // Продолжение операции

                        // Вместо сложной операции, используя std::complex<double>
                        std::complex<double> exp_factor_rx = exp(j) * exp(2 * M_PI * rx / 0.1);
                        std::complex<double> exp_factor_tx = exp(j) * exp(2 * M_PI * tx / 0.1);
                        std::complex<double> channelCoefficients_n = temp2(0, 0) * exp_factor_rx * exp_factor_tx;

                        channelCoefficients_u_s_n(s + u + pair, n + 4) += channelCoefficients_n;

                    }
                    else if (n < 2) {
                        Eigen::Vector2d F_tx = transmitter.FieldPattern(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), ksi_tx);
                        Eigen::Vector2d F_rx = receiver.FieldPattern(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), ksi_rx);

                        //Eigen::Vector2d F_tx = transmitter.transformationFromLCStoGCS(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F1_tx);
                        //Eigen::Vector2d F_rx = receiver.transformationFromLCStoGCS(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F1_rx);

                        Eigen::Vector3d sphericalUnitVector_tx(sin(thetaZOD_n_m(n, m)) * cos(phiAOD_n_m(n, m)),
                            sin(thetaZOD_n_m(n, m)) * sin(phiAOD_n_m(n, m)),
                            cos(thetaZOD_n_m(n, m)));
                        Eigen::Vector3d sphericalUnitVector_rx(sin(thetaZOA_n_m(n, m)) * cos(phiAOA_n_m(n, m)),
                            sin(thetaZOA_n_m(n, m)) * sin(phiAOA_n_m(n, m)),
                            cos(thetaZOA_n_m(n, m)));


                        Eigen::Matrix2cd XPR_and_InitialRandomPhases;
                        XPR_and_InitialRandomPhases <<
                            exp(j * initialPhases(n, m * 4)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 1)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 2)),
                            exp(j * initialPhases(n, m * 4 + 3));


                        double tx = sphericalUnitVector_tx.transpose() * transmitter.generateAntennaElements().row(s).transpose();
                        double rx = sphericalUnitVector_rx.transpose() * receiver.generateAntennaElements().row(u).transpose();



                        // Разделение сложной операции на более простые
                        auto temp1 = F_rx.transpose() * XPR_and_InitialRandomPhases; // Промежуточный результат
                        auto temp2 = temp1 * F_tx; // Продолжение операции


                        // Вместо сложной операции, используя std::complex<double>
                        std::complex<double> exp_factor_rx = exp(j) * exp(2 * M_PI * rx / 0.1);
                        std::complex<double> exp_factor_tx = exp(j) * exp(2 * M_PI * tx / 0.1);
                        std::complex<double> channelCoefficients_n = temp2(0, 0) * exp_factor_rx * exp_factor_tx;


                        if (m < 8 || m >17) {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n) += 0.05 * channelCoefficients_n;
                        }
                        else if ((m > 8 && m < 12) || (m > 15 && m < 18)) {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n + 1) += 0.05 * channelCoefficients_n;
                        }
                        else {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n + 2) += 0.05 * channelCoefficients_n;
                        }

                    }
                }
            }
            pair = pair + 15;
        }
        channelCoefficients_u_s_n.col(n) *= sqrt(clusterPowers[n]);
    }

    return channelCoefficients_u_s_n;
}

//ИХ для ЛОС
Eigen::MatrixXcd generateLOSChannelCoefficients(const UserTerminal& transmitter, const UserTerminal& receiver,
    std::vector<double>& clusterPowers, Eigen::MatrixXd& phiAOD_n_m, Eigen::MatrixXd& phiAOA_n_m,
    Eigen::MatrixXd& thetaZOD_n_m, Eigen::MatrixXd& thetaZOA_n_m, Eigen::MatrixXd& XRP, Eigen::MatrixXd& initialPhases, double riceanK) {

    riceanK = pow(10, riceanK / 10);
    phiAOD_n_m = phiAOD_n_m * M_PI / 180;
    phiAOA_n_m = phiAOA_n_m * M_PI / 180;
    thetaZOD_n_m = thetaZOD_n_m * M_PI / 180;
    thetaZOA_n_m = thetaZOA_n_m * M_PI / 180;

    std::complex<double> j(0.0, 1.0);
    Eigen::MatrixXcd channelCoefficients_u_s_n(256, clusterPowers.size() + 4);//16*16 комбинаций t-u
    channelCoefficients_u_s_n.setZero();

    // пары U-S
    for (int n = 0; n < clusterPowers.size(); ++n) {
        double ksi_rx;
        double ksi_tx;
        int pair = 0;
        for (int u = 0; u < 16; ++u) {
            if (u < 8) { ksi_rx = 45; }
            else { ksi_rx = -45; }
            for (int s = 0; s < 16; ++s) {
                if (s < 8) { ksi_tx = 45; }
                else { ksi_tx = -45; }
                for (int m = 0; m < 20; ++m) {
                    if (n > 2) {
                        Eigen::Vector2d F_tx = transmitter.FieldPattern(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), ksi_tx);
                        Eigen::Vector2d F_rx = receiver.FieldPattern(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), ksi_rx);

                        //Eigen::Vector2d F_tx = transmitter.transformationFromLCStoGCS(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F1_tx);
                        //Eigen::Vector2d F_rx = receiver.transformationFromLCStoGCS(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F1_rx);

                        Eigen::Vector3d sphericalUnitVector_tx(sin(thetaZOD_n_m(n, m)) * cos(phiAOD_n_m(n, m)),
                            sin(thetaZOD_n_m(n, m)) * sin(phiAOD_n_m(n, m)),
                            cos(thetaZOD_n_m(n, m)));
                        Eigen::Vector3d sphericalUnitVector_rx(sin(thetaZOA_n_m(n, m)) * cos(phiAOA_n_m(n, m)),
                            sin(thetaZOA_n_m(n, m)) * sin(phiAOA_n_m(n, m)),
                            cos(thetaZOA_n_m(n, m)));

                        Eigen::Matrix2cd XPR_and_InitialRandomPhases;
                        XPR_and_InitialRandomPhases <<
                            exp(j * initialPhases(n, m * 4)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 1)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 2)),
                            exp(j * initialPhases(n, m * 4 + 3));

                        double tx = sphericalUnitVector_tx.transpose() * transmitter.generateAntennaElements().row(s).transpose();
                        double rx = sphericalUnitVector_rx.transpose() * receiver.generateAntennaElements().row(u).transpose();

                        // Разделение сложной операции на более простые
                        auto temp1 = F_rx.transpose() * XPR_and_InitialRandomPhases; // Промежуточный результат
                        auto temp2 = temp1 * F_tx; // Продолжение операции

                        // Вместо сложной операции, используя std::complex<double>
                        std::complex<double> exp_factor_rx = exp(j) * exp(2 * M_PI * rx / 0.1);
                        std::complex<double> exp_factor_tx = exp(j) * exp(2 * M_PI * tx / 0.1);
                        std::complex<double> channelCoefficients_n = temp2(0, 0) * exp_factor_rx * exp_factor_tx;

                        channelCoefficients_u_s_n(s + u + pair, n + 4) += channelCoefficients_n;
                    }
                    else if (n > 0 && n < 3) {
                        Eigen::Vector2d F_tx = transmitter.FieldPattern(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), ksi_tx);
                        Eigen::Vector2d F_rx = receiver.FieldPattern(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), ksi_rx);

                        //Eigen::Vector2d F_tx = transmitter.transformationFromLCStoGCS(thetaZOD_n_m(n, m), phiAOD_n_m(n, m), transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F1_tx);
                        //Eigen::Vector2d F_rx = receiver.transformationFromLCStoGCS(thetaZOA_n_m(n, m), phiAOA_n_m(n, m), receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F1_rx);

                        Eigen::Vector3d sphericalUnitVector_tx(sin(thetaZOD_n_m(n, m)) * cos(phiAOD_n_m(n, m)),
                            sin(thetaZOD_n_m(n, m)) * sin(phiAOD_n_m(n, m)),
                            cos(thetaZOD_n_m(n, m)));
                        Eigen::Vector3d sphericalUnitVector_rx(sin(thetaZOA_n_m(n, m)) * cos(phiAOA_n_m(n, m)),
                            sin(thetaZOA_n_m(n, m)) * sin(phiAOA_n_m(n, m)),
                            cos(thetaZOA_n_m(n, m)));

                        Eigen::Matrix2cd XPR_and_InitialRandomPhases;
                        XPR_and_InitialRandomPhases <<
                            exp(j * initialPhases(n, m * 4)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 1)),
                            sqrt(1 / XRP(n, m))* exp(j * initialPhases(n, m * 4 + 2)),
                            exp(j * initialPhases(n, m * 4 + 3));

                        double tx = sphericalUnitVector_tx.transpose() * transmitter.generateAntennaElements().row(s).transpose();
                        double rx = sphericalUnitVector_rx.transpose() * receiver.generateAntennaElements().row(u).transpose();

                        // Разделение сложной операции на более простые
                        auto temp1 = F_rx.transpose() * XPR_and_InitialRandomPhases; // Промежуточный результат
                        auto temp2 = temp1 * F_tx; // Продолжение операции

                        // Вместо сложной операции, используя std::complex<double>
                        std::complex<double> exp_factor_rx = exp(j) * exp(2 * M_PI * rx / 0.1);
                        std::complex<double> exp_factor_tx = exp(j) * exp(2 * M_PI * tx / 0.1);
                        std::complex<double> channelCoefficients_n = temp2(0, 0) * exp_factor_rx * exp_factor_tx;

                        if (m < 8 || m >17) {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n - 2) += 0.05 * channelCoefficients_n;
                        }
                        else if ((m > 8 && m < 12) || (m > 15 && m < 18)) {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n + -1) += 0.05 * channelCoefficients_n;
                        }
                        else {
                            channelCoefficients_u_s_n(s + u + pair, 3 * n) += 0.05 * channelCoefficients_n;
                        }

                    }

                }

                if (n == 0) {
                    Eigen::Vector2d F_tx = transmitter.FieldPattern(thetaZOD_n_m(n, 0), phiAOD_n_m(n, 0), 0.0);
                    Eigen::Vector2d F_rx = transmitter.FieldPattern(thetaZOA_n_m(n, 0), phiAOA_n_m(n, 0), 0.0);

                    //Eigen::Vector2d F_tx = transmitter.transformationFromLCStoGCS(thetaZOD_n_m(n, 0), phiAOD_n_m(n, 0), transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F1_tx);
                    //Eigen::Vector2d F_rx = receiver.transformationFromLCStoGCS(thetaZOA_n_m(n, 0), phiAOA_n_m(n, 0), receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F1_rx);

                    Eigen::Vector3d sphericalUnitVector_tx(sin(thetaZOD_n_m(n, 0)) * cos(phiAOD_n_m(n, 0)),
                        sin(thetaZOD_n_m(n, 0)) * sin(phiAOD_n_m(n, 0)),
                        cos(thetaZOD_n_m(n, 0)));
                    Eigen::Vector3d sphericalUnitVector_rx(sin(thetaZOA_n_m(n, 0)) * cos(phiAOA_n_m(n, 0)),
                        sin(thetaZOA_n_m(n, 0)) * sin(phiAOA_n_m(n, 0)),
                        cos(thetaZOA_n_m(n, 0)));


                    Eigen::Matrix2cd XPR_and_InitialRandomPhases;
                    XPR_and_InitialRandomPhases << 1.0, 0.0, 0.0, -1.0;


                    double tx = sphericalUnitVector_tx.transpose() * transmitter.generateAntennaElements().row(s).transpose();
                    double rx = sphericalUnitVector_rx.transpose() * receiver.generateAntennaElements().row(u).transpose();

                    // Разделение сложной операции на более простые
                    auto temp1 = F_rx.transpose() * XPR_and_InitialRandomPhases; // Промежуточный результат
                    auto temp2 = temp1 * F_tx; // Продолжение операции


                    // Вместо сложной операции, используя std::complex<double>
                    std::complex<double> exp_factor_rx = exp(j) * exp(2 * M_PI * rx / 0.1);
                    std::complex<double> exp_factor_tx = exp(j) * exp(2 * M_PI * tx / 0.1);
                    std::complex<double> exp_d3D = exp(j) * exp(-2 * M_PI * calculateDistance(transmitter, receiver) / ( 0.1));

                    std::complex<double> channelCoefficients_n = temp2(0, 0) * exp_factor_rx * exp_factor_tx * exp_d3D;

                    channelCoefficients_u_s_n(s + u + pair, n) += channelCoefficients_n;
                }
            }
            pair = pair + 15;
        }
        if (n != 0) {
            channelCoefficients_u_s_n.col(n) *= sqrt(clusterPowers[n]);
            channelCoefficients_u_s_n.col(n) *= sqrt(1 / (riceanK + 1));
        }
        else {
            channelCoefficients_u_s_n.col(n) *= sqrt(riceanK / (riceanK + 1));
        }
    }

    return channelCoefficients_u_s_n;
}

//____________________________________________Основная_Программа___________________________________________//
int main() {
    // Параметры комнаты 
    double roomLength = 120.0; // Длина комнаты
    double roomWidth = 50.0;    // Ширина комнаты 
    double roomHeight = 3.0;    // Высота комнаты 
    double wavelength = 0.3;     // Длина волны 

    // Создаем пользователей 
    std::vector<UserTerminal> users;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> xDist(0.0, roomLength);
    std::uniform_real_distribution<> yDist(0.0, roomWidth);
    double userHeight = 1.0; // Высота пользователей 

    for (int i = 1; i <= 12; ++i) {
        //double bearing = (rand() % 90) * M_PI / 180.0; // Угол поворота 
        double bearing = (std::rand() % 360) * M_PI / 180.0;
        double downtilt = (rand() % 90) * M_PI / 180.0; // Угол наклона
        double slant = (rand() % 90) * M_PI / 180.0; // Угол наклона

        UserTerminal newUT(i, 0, 0, userHeight, bearing, downtilt, slant);
        bool isValidPosition = false;

        while (!isValidPosition) {
            newUT.x = xDist(gen);
            newUT.y = yDist(gen);

            // Проверяем расстояние до всех существующих пользователей
            isValidPosition = true; // Предполагаем, что позиция валидна
            for (const auto& user : users) {
                if (calculateDistance(newUT, user) < 1.0) {
                    isValidPosition = false; // Позиция невалидна
                    break; // Выходим из цикла
                }
            }
            users.push_back(newUT);
        }
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

        std::cout << "// UT " << transmitter.id << " - UT " << receiver.id << " // " << std::endl << std::endl;

        // Расчёт лос угов и соовтетствующих значений ДН полей приёмника и передатчика
        double losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA;
        transmitter.calculateLOSAngles(transmitter, receiver, losPhiAOD, losThetaZOD, losPhiAOA, losThetaZOA);

        std::cout << "LOS ZOD, AOD: (" << losThetaZOD << ", " << losPhiAOD << "),\n"
            << "LOS ZOA, AOA: (" << losThetaZOA << ", " << losPhiAOA << ")\n\n";

        Eigen::Vector2d F_tx = transmitter.FieldPattern(losThetaZOD, losPhiAOD, 45);
        std::cout << "{F_tx_theta,F_tx_pfi} : " << F_tx[0] << " ; " << F_tx[1] << std::endl;
        /*
        Eigen::Vector2d txAntennaPattern = transmitter.transformationFromLCStoGCS(losThetaZOD, losPhiAOD, transmitter.bearingAngle, transmitter.downtiltAngle, transmitter.slantAngle, F_tx);
        std::cout << "Bearing Angle for transmitter  = " << transmitter.downtiltAngle << " rad" << std::endl;
        std::cout << "Transformation from LCS to GCS F_tx_Theta, F_tx_Pfi  : { " << txAntennaPattern[0] << " ; " << txAntennaPattern[1] << " }" << std::endl << std::endl;
        */

        Eigen::Vector2d F_rx = receiver.FieldPattern(losThetaZOA, losPhiAOA, 45);
        std::cout << "{F_rx_theta,F_rx_pfi} : " << F_rx[0] << " ; " << F_rx[1] << std::endl;
        /*
        Eigen::Vector2d rxAntennaPattern = receiver.transformationFromLCStoGCS(losThetaZOD, losPhiAOD, receiver.bearingAngle, receiver.downtiltAngle, receiver.slantAngle, F_rx);
        std::cout << "Bearing Angle for receiver  = " << receiver.downtiltAngle << " rad" << std::endl;
        std::cout << "Transformation from LCS to GCS F_rx_Theta, F_rx_Pfi  : { " << rxAntennaPattern[0] << " ; " << rxAntennaPattern[1] << " }" << std::endl << std::endl;
        */

        // Элементы антенны
        Eigen::MatrixXd d_tx = transmitter.generateAntennaElements();
        std::cout << "d_tx :\n" << d_tx << std::endl << std::endl;

        //_______________________________________STEP_2________________________________________//

        // Вычисление расстояния между передатчиком и приемником
        double d = calculateDistance(transmitter, receiver);
        std::cout << "Distance between transmitter and receiver is " << d << std::endl;

        double P_LOS;
        bool los;
        if (d <= 5)
        {
            P_LOS = 1;
            los = true;
            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
        }
        else if (d > 5 && d <= 49)
        {
            P_LOS = exp(-(d - 5) / 70.8);
            los = true;
            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
        }
        else
        {
            P_LOS = exp(-(d - 49) / 211.7) * 0.54;
            if (P_LOS < 0.5) {
                los = false;
                std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - NLOS" << std::endl;
            }
            else {
                los = true;
                std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
            }
        }

        //_________________________________________________STEP_3__________________________________________//
        // Вычисление Pathloss_LOS для каждой пары
        double path_loss_LOS;
        double nu = 3.0; // частота в ГГц
        path_loss_LOS = 32.4 + 17.3 * log10(d) + 20 * log10(nu);

        std::cout << "PathLoss in LOS between users " << transmitter.id << " and " << receiver.id << " = " << path_loss_LOS << std::endl << std::endl;

        //_________________________________________________STEP_4__________________________________________//
        std::cout << "\nLSP for LOS for User " << transmitter.id << " and User " << receiver.id << " :\n\n";
        LargeScaleParameters lsp(los, nu); //3ГГц
        lsp.showParameters();

        //______________STEP_5_______________//
        // Генерация задержек кластеров
        std::vector<double> clusterDelays = generateClusterDelays(los, lsp.delaySpread,  lsp.riceanK); // Передаем delaySpread из LSP
        
        //_____________STEP_6_______________//
        // Генерация мощностей кластеров
        std::vector<double> clusterPowers = generateClusterPowers(los, clusterDelays, lsp.delaySpread);

        //оставляю лишь нужные задержки  и сортирую их с мощностями 
        for (int i = indicesToDelete.size() - 1; i >= 0; --i) {
            clusterDelays.erase(clusterDelays.begin() + indicesToDelete[i]);
        }
        indicesToDelete.clear();
        std::pair<std::vector<double>, std::vector<double>> sortClusterPowers_ClusterDelays = sort_with_indices(clusterPowers, clusterDelays);
        clusterPowers = sortClusterPowers_ClusterDelays.first;
        clusterDelays = sortClusterPowers_ClusterDelays.second;

        std::vector<double> clusterPowersWithScalingFactors = clusterPowers;
        if (los) {

            for (size_t n = 0; n < clusterPowersWithScalingFactors.size(); ++n) {
                clusterPowersWithScalingFactors[n] *= (1 / (pow(10, lsp.riceanK) + 1))  ;
            }
            clusterPowersWithScalingFactors[0] += (pow(10, lsp.riceanK) / (pow(10, lsp.riceanK) + 1));
        }


        // Выводы значений 
        std::cout << "Cluster delays for User " << transmitter.id << " and User " << receiver.id << ":\n\n";
        int j = 1;
        for (const auto& delay : clusterDelays) {
            std::cout << "-delay for power #" << j << " : " << delay << "\n";
            j++;
        }
        std::cout << std::endl;

        std::cout << "Cluster powers for User " << transmitter.id << " and User " << receiver.id << ": ";
        for (const auto& power : clusterPowers) {
            std::cout << power << " ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Cluster powers with scaling factors for User " << transmitter.id << " and User " << receiver.id << ": ";
        for (const auto& powerWithScalingFactors : clusterPowersWithScalingFactors) {
            std::cout << powerWithScalingFactors << " ";
        }
        std::cout << std::endl << std::endl;

        //_____________STEP_7_______________//

        //Generate arrival angles and departure angles for both azimuth and elevation
        Eigen::MatrixXd PhiAOD = generatePhiAOD(los, clusterPowersWithScalingFactors, lsp.azimuthSpreadDeparture, lsp.riceanK, losPhiAOD);
        Eigen::MatrixXd PhiAOA = generatePhiAOA(los, clusterPowersWithScalingFactors, lsp.azimuthSpreadArrival, lsp.riceanK, losPhiAOA);
        Eigen::MatrixXd ThetaZOD = generateThetaZOD(los, clusterPowersWithScalingFactors, lsp.zenithSpreadDeparture, lsp.riceanK, losThetaZOD);
        Eigen::MatrixXd ThetaZOA = generateThetaZOA(los, clusterPowersWithScalingFactors, lsp.zenithSpreadArrival, lsp.riceanK, losThetaZOA);

        if (los) {
            std::cout << "PhiAOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << PhiAOD(0, 0) << std::endl << PhiAOD.bottomRows(PhiAOD.rows() - 1) << std::endl << std::endl;
            std::cout << "PhiAOA for UT receiver " << receiver.id << ": \n";
            std::cout << PhiAOA(0, 0) << std::endl << PhiAOA.bottomRows(PhiAOA.rows() - 1) << std::endl << std::endl;
            std::cout << "ThetaZOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << ThetaZOD(0, 0) << std::endl << ThetaZOD.bottomRows(ThetaZOD.rows() - 1) << std::endl << std::endl;
            std::cout << "ThetaZOA for UT receiver " << receiver.id << ": \n";
            std::cout << ThetaZOA(0, 0) << std::endl << ThetaZOA.bottomRows(ThetaZOA.rows() - 1) << std::endl << std::endl;
        }
        else {
            std::cout << "PhiAOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << PhiAOD << std::endl << std::endl;
            std::cout << "PhiAOA for UT receiver " << receiver.id << ": \n";
            std::cout << PhiAOA << std::endl << std::endl;
            std::cout << "ThetaZOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << ThetaZOD << std::endl << std::endl;
            std::cout << "ThetaZOA for UT receiver " << receiver.id << ": \n";
            std::cout << ThetaZOA << std::endl << std::endl;
        }



        std::cout << "Coupling of rays within a cluster for both azimuth and elevation " << std::endl;
        randomCouplingRays(PhiAOD, PhiAOA, ThetaZOD, ThetaZOA, los);

        if (los) {
            std::cout << "PhiAOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << PhiAOD(0, 0) << std::endl << PhiAOD.bottomRows(PhiAOD.rows() - 1) << std::endl << std::endl;
            std::cout << "PhiAOA for UT receiver " << receiver.id << ": \n";
            std::cout << PhiAOA(0, 0) << std::endl << PhiAOA.bottomRows(PhiAOA.rows() - 1) << std::endl << std::endl;
            std::cout << "ThetaZOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << ThetaZOD(0, 0) << std::endl << ThetaZOD.bottomRows(ThetaZOD.rows() - 1) << std::endl << std::endl;
            std::cout << "ThetaZOA for UT receiver " << receiver.id << ": \n";
            std::cout << ThetaZOA(0, 0) << std::endl << ThetaZOA.bottomRows(ThetaZOA.rows() - 1) << std::endl << std::endl;
        }
        else {
            std::cout << "PhiAOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << PhiAOD << std::endl << std::endl;
            std::cout << "PhiAOA for UT receiver " << receiver.id << ": \n";
            std::cout << PhiAOA << std::endl << std::endl;
            std::cout << "ThetaZOD for UT transmitter " << transmitter.id << ": \n";
            std::cout << ThetaZOD << std::endl << std::endl;
            std::cout << "ThetaZOA for UT receiver " << receiver.id << ": \n";
            std::cout << ThetaZOA << std::endl << std::endl;
        }

        //___________A1___A2___________//
        std::vector<double> AS = calculateAngularSpreadandMeanAngles(los, clusterPowers, PhiAOD, PhiAOA, ThetaZOD, ThetaZOA);
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






        //_____________STEP_9_______________//

        //XPR
        Eigen::MatrixXd XPR = generateXPR(los, clusterPowers);
        std::cout << "Generate the cross polarization power ratios  K_n_m: \n";
        std::cout << XPR << std::endl << std::endl;

        //_____________STEP_10_______________//
        Eigen::MatrixXd initialPhases = generateInitialRandomPhases(clusterPowers);

        std::cout << "Initial random phases for each ray in each cluster:" << std::endl;
        for (int n = 0; n < clusterPowers.size(); ++n) {
            std::cout << "Cluster " << n + 1 << ": \n";
            for (int m = 0; m < 20; ++m) {
                std::cout << "Ray " << m + 1 << " (Theta - Theta, Theta - Phi, Phi - Theta, Phi - Phi): | "
                    << initialPhases(n, m * 4) << ", "
                    << initialPhases(n, m * 4 + 1) << ", "
                    << initialPhases(n, m * 4 + 2) << ", "
                    << initialPhases(n, m * 4 + 3) << " |\n";
            }
            std::cout << std::endl;
        }

        if (!los) {
            Eigen::MatrixXcd channelСoefficients = generateNLOSChannelCoefficients(transmitter, receiver, clusterPowers, PhiAOD, PhiAOA, ThetaZOD, ThetaZOA, XPR, initialPhases);
            Eigen::MatrixXd modulusMatrix = channelСoefficients.array().abs();
            std::cout << "Module channelCoefficients:\n" << modulusMatrix << std::endl;
        }
        else {
            Eigen::MatrixXcd channelСoefficients = generateLOSChannelCoefficients(transmitter, receiver, clusterPowers, PhiAOD, PhiAOA, ThetaZOD, ThetaZOA, XPR, initialPhases, lsp.riceanK);
            Eigen::MatrixXd modulusMatrix = channelСoefficients.array().abs();
            std::cout << "Module channelCoefficients:\n" << modulusMatrix << std::endl;
        }
    }
    return 0;
}