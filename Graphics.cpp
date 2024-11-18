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

#include <string>
#include <fstream>
#include <filesystem>

using namespace Eigen;

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
            MatrixXd C(7, 7);

            //Если LOS :
            //Затемнение SF
            StandardDeviationSF = 3;
            MeanSF = 0;
            /*std::normal_distribution<> loSshadowFadingDist(MeanSF, StandardDeviationSF);
            shadowFading = loSshadowFadingDist(gen);*/

            //К-фактор (K)
            StandardDeviationK = 4;
            MeanK = 7;
            /*std::normal_distribution<> riceanKDist0(MeanK, StandardDeviationK);
            riceanK = (riceanKDist0(gen));*/

            // разброс задержки (DS)
            StandardDeviationDS = 0.18;
            MeanDS = (-0.01 * log10(1 + fc) - 7.692);
            /*std::normal_distribution<> loSdelaySpreadDist(MeanDS, StandardDeviationDS);
            delaySpread = loSdelaySpreadDist(gen);*/

            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.18;
            MeanASD = 1.6;
            /*std::normal_distribution<> loSAzimuthSpreadDepartureDist(MeanASD, StandardDeviationASD);
            azimuthSpreadDeparture = std::min((loSAzimuthSpreadDepartureDist(gen)), log10(104.0));*/

            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.119);
            MeanASA = (-0.19 * log10(1 + fc) + 1.781);
            /*std::normal_distribution<> loSAzimuthSpreadArrivalDist(MeanASA, StandardDeviationASA);
            azimuthSpreadArrival = std::min((loSAzimuthSpreadDepartureDist(gen)), log10(104.0));*/


            //For frequencies below 6 GHz, use fc = 6 when determining the values of the frequency-dependent
            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.04 * log10(1 + fc) + 0.264);
            MeanZSA = (-0.26 * log10(1 + fc) + 1.44);
            /*std::normal_distribution<> loSZenithSpreadArrivalDist(MeanZSA, StandardDeviationZSA);
            zenithSpreadArrival = std::min(loSZenithSpreadArrivalDist(gen), log10(52.0));*/

            //Зенитный угол Разброс вылета (ZSD)

            if (fc < 6.0) {
                StandardDeviationZSD = (0.13 * log10(1 + 6) + 0.30);
                MeanZSD = (-1.43 * log10(1 + 6) + 2.228);
                /*std::normal_distribution<> loSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
                zenithSpreadDeparture = std::min(loSZenithSpreadDepartureDist(gen), log10(52.0));*/
            }
            else {
                StandardDeviationZSD = (0.13 * log10(1 + fc) + 0.30);
                MeanZSD = (-1.43 * log10(1 + fc) + 2.228);
                /*std::normal_distribution<> loSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
                zenithSpreadDeparture = std::min(loSZenithSpreadDepartureDist(gen), log10(52.0));*/
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

            std::cout << value << "\n";
            /*
            losStandardDeviationSF = pow(10, losStandardDeviationSF/10);
            losStandardDeviationK = pow(10, losStandardDeviationK / 10);
            losStandardDeviationDS = pow(10, losStandardDeviationDS );
            losStandardDeviationASD = pow(10, losStandardDeviationASD );
            losStandardDeviationASA = pow(10, losStandardDeviationASA);
            losStandardDeviationZSD = pow(10, losStandardDeviationZSD );
            losStandardDeviationZSA = pow(10, losStandardDeviationZSA);
            */

            //SF K DS ASD ASA ZSD ZSA
            C << StandardDeviationSF * StandardDeviationSF, 0.5 * StandardDeviationSF * StandardDeviationK, -0.8 * StandardDeviationSF * StandardDeviationDS, -0.4 * StandardDeviationSF * StandardDeviationASD, -0.5 * StandardDeviationSF * StandardDeviationASA, 0.2 * StandardDeviationSF * StandardDeviationZSD, 0.3 * StandardDeviationSF * StandardDeviationZSA,
                0.5 * StandardDeviationK * StandardDeviationSF, StandardDeviationK* StandardDeviationK, -0.5 * StandardDeviationK * StandardDeviationDS, 0.0 * StandardDeviationK * StandardDeviationASD, 0.0 * StandardDeviationK * StandardDeviationASA, 0.0 * StandardDeviationK * StandardDeviationZSD, 0.1 * StandardDeviationK * StandardDeviationZSA,
                -0.8 * StandardDeviationDS * StandardDeviationSF, -0.5 * StandardDeviationDS * StandardDeviationK, StandardDeviationDS* StandardDeviationDS, 0.6 * StandardDeviationDS * StandardDeviationASD, 0.8 * StandardDeviationDS * StandardDeviationASA, 0.1 * StandardDeviationDS * StandardDeviationZSD, 0.2 * StandardDeviationDS * StandardDeviationZSA,
                -0.4 * StandardDeviationASD * StandardDeviationSF, 0.0 * StandardDeviationASD * StandardDeviationK, 0.6 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.4 * StandardDeviationASD * StandardDeviationASA, 0.5 * StandardDeviationASD * StandardDeviationZSD, 0.0 * StandardDeviationASD * StandardDeviationZSA,
                -0.5 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationK, 0.8 * StandardDeviationASA * StandardDeviationDS, 0.4 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, 0.0 * StandardDeviationASA * StandardDeviationZSD, 0.5 * StandardDeviationASA * StandardDeviationZSA,
                0.2 * StandardDeviationZSD * StandardDeviationSF, 0.0 * StandardDeviationZSD * StandardDeviationK, 0.1 * StandardDeviationZSD * StandardDeviationDS, 0.5 * StandardDeviationZSD * StandardDeviationASD, 0.0 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.0 * StandardDeviationZSD * StandardDeviationZSA,
                0.3 * StandardDeviationZSA * StandardDeviationSF, 0.1 * StandardDeviationZSA * StandardDeviationK, 0.2 * StandardDeviationZSA * StandardDeviationDS, 0.0 * StandardDeviationZSA * StandardDeviationASD, 0.5 * StandardDeviationZSA * StandardDeviationASA, 0.0 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;


            std::cout << C << std::endl << std::endl;


            // Вычисление и проверка главных миноров
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


            MatrixXd L;
            L = C.llt().matrixL();
            std::cout << L << "\n";

            // Проверка на положительную определенность
            if (C.llt().info() != Success) {
                std::cerr << "Matrix C is not positive definite!" << std::endl;
                return;
            }

            // Генерация стандартных нормальных случайных величин
            std::mt19937 gen;
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
            MatrixXd C(6, 6);

            //Если NLOS:
            //Затенение (SF)
            StandardDeviationSF = 8.03;
            MeanSF = 0;
            /*std::normal_distribution<> nLoSshadowFadingDist(MeanSF, StandardDeviationSF);
            shadowFading = nLoSshadowFadingDist(gen);*/

            // разброс задержки (DS)
            StandardDeviationDS = (0.10 * log10(1 + fc) + 0.055);
            MeanDS = (-0.28 * log10(1 + fc) - 7.173);
            /*std::normal_distribution<> nLoSdelaySpreadDist(MeanDS, StandardDeviationDS);
            delaySpread = nLoSdelaySpreadDist(gen);*/ // Задержка распространения НЛОС

            //Азимутальный угол Разброс вылета (ASD)
            StandardDeviationASD = 0.25;
            MeanASD = 1.62;
            /*std::normal_distribution<> nLoSAzimuthSpreadDepartureDist(MeanASD, StandardDeviationASD);
            azimuthSpreadDeparture = std::min((nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));*/

            //Азимутальный угол Разброс прихода (ASA)
            StandardDeviationASA = (0.12 * log10(1 + fc) + 0.059);
            MeanASA = (-0.11 * log10(1 + fc) + 1.863);
            /*std::normal_distribution<> nLoSAzimuthSpreadArrivaDist(MeanASA, StandardDeviationASA);
            azimuthSpreadArrival = std::min((nLoSAzimuthSpreadDepartureDist(gen)), log10(104.0));*/


            //Зенитный угол Разброс прихода (ZSA)
            StandardDeviationZSA = (-0.09 * log10(1 + fc) + 0.746);
            MeanZSA = (-0.15 * log10(1 + fc) + 1.387);
            /*std::normal_distribution<> nLoSZenithSpreadArrivalDist(MeanZSA, StandardDeviationZSA);
            zenithSpreadArrival = std::min(nLoSZenithSpreadArrivalDist(gen), log10(52.0));*/

            //Зенитный угол Разброс вылета (ZSD)
            StandardDeviationZSD = 0.36;
            MeanZSD = 1.08;

            /*std::normal_distribution<> nLoSZenithSpreadDepartureDist(MeanZSD, StandardDeviationZSD);
            zenithSpreadDeparture = std::min(nLoSZenithSpreadDepartureDist(gen), log10(52.0));*/

            value << shadowFading, delaySpread, azimuthSpreadDeparture, azimuthSpreadArrival, zenithSpreadDeparture, zenithSpreadArrival;
            means << MeanSF, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

            std::cout << value << "\n";
            /*
            losStandardDeviationSF = pow(10, losStandardDeviationSF/10);
            losStandardDeviationK = pow(10, losStandardDeviationK / 10);
            losStandardDeviationDS = pow(10, losStandardDeviationDS );
            losStandardDeviationASD = pow(10, losStandardDeviationASD );
            losStandardDeviationASA = pow(10, losStandardDeviationASA);
            losStandardDeviationZSD = pow(10, losStandardDeviationZSD );
            losStandardDeviationZSA = pow(10, losStandardDeviationZSA);
            */

            //SF  DS ASD ASA ZSD ZSA
            C << StandardDeviationSF * StandardDeviationSF, -0.5 * StandardDeviationSF * StandardDeviationDS, 0.0 * StandardDeviationSF * StandardDeviationASD, -0.4 * StandardDeviationSF * StandardDeviationASA, 0.0 * StandardDeviationSF * StandardDeviationZSD, 0.0 * StandardDeviationSF * StandardDeviationZSA,
                -0.5 * StandardDeviationDS * StandardDeviationSF, StandardDeviationDS* StandardDeviationDS, 0.4 * StandardDeviationDS * StandardDeviationASD, 0.0 * StandardDeviationDS * StandardDeviationASA, -0.27 * StandardDeviationDS * StandardDeviationZSD, -0.06 * StandardDeviationDS * StandardDeviationZSA,
                0.0 * StandardDeviationASD * StandardDeviationSF, 0.4 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.0 * StandardDeviationASD * StandardDeviationASA, 0.35 * StandardDeviationASD * StandardDeviationZSD, 0.23 * StandardDeviationASD * StandardDeviationZSA,
                -0.4 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationDS, 0.0 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, -0.08 * StandardDeviationASA * StandardDeviationZSD, 0.43 * StandardDeviationASA * StandardDeviationZSA,
                0.0 * StandardDeviationZSD * StandardDeviationSF, -0.27 * StandardDeviationZSD * StandardDeviationDS, 0.35 * StandardDeviationZSD * StandardDeviationASD, -0.08 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.42 * StandardDeviationZSD * StandardDeviationZSA,
                0.0 * StandardDeviationZSA * StandardDeviationSF, -0.06 * StandardDeviationZSA * StandardDeviationDS, 0.23 * StandardDeviationZSA * StandardDeviationASD, 0.43 * StandardDeviationZSA * StandardDeviationASA, 0.42 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

            std::cout << C << std::endl << std::endl;

            // Вычисление и проверка главных миноров
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

            MatrixXd L;
            L = C.llt().matrixL();
            std::cout << L << "\n";

            // Проверка на положительную определенность
            if (C.llt().info() != Success) {
                std::cerr << "Matrix C is not positive definite!" << std::endl;
                return;
            }

            // Генерация стандартных нормальных случайных величин
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

//___________________________________________________Класс_UT______________________________________________//
class UserTerminal {
public:
    int id;
    double x, y, z; // Координаты пользователя 
    double wavelength; // Длина волны 
    int numElementsX = 4; // 4 элемента антенны по X
    int numElementsY = 2; // 2 элемента антенны по Y 
    double bearingAngle; // Угол поворота 
    double downtiltAngle; // Угол наклона 
    double slantAngle; // Угол наклона

    // Конструктор класса

    UserTerminal(int id, double x, double y, double z, double lambda, double bearing, double downtilt, double slant)
        : id(id), x(x), y(y), z(z), wavelength(lambda), bearingAngle(bearing), downtiltAngle(downtilt), slantAngle(slant) {
    }

    // Методы для вычисления ДН полей
    Vector2d FieldPattern(double thetaAngle, double phiAngle) const {
        Vector2d fieldPattern;
        thetaAngle = thetaAngle * 180 / M_PI;
        phiAngle = phiAngle * 180 / M_PI;

        double ksi = 45.0 * M_PI / 180;

        //Считалось ранее по модели 1 
        //double f2_Theta = sqrt(std::min(12 * (phiAngle / 65) * (phiAngle / 65), 30.0));
        //double f2_Phi = sqrt(std::min(12 * ((thetaAngle - 90) / 65) * ((thetaAngle - 90) / 65), 30.0));
        //double cos_Pci = (cos(ksi) * sin(thetaAngle) + sin(ksi) * sin(phiAngle) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));
        //double sin_Pci = (sin(ksi) * cos(thetaAngle)) / (pow((1 - (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle)) * (cos(ksi) * cos(thetaAngle) - sin(ksi) * sin(phiAngle) * cos(thetaAngle))), 0.5));
        //double f1_Theta = sin_Pci * f2_Theta + cos_Pci * f2_Phi;
        //double f1_Phi = cos_Pci * f2_Theta - sin_Pci * f2_Phi;


        double AverticalPowerPattern = (-1) * std::min(12 * (phiAngle / 65) * (phiAngle / 65), 30.0);
        double AhorizontalPowerPattern = (-1) * std::min(12 * ((thetaAngle - 90) / 65) * ((thetaAngle - 90) / 65), 30.0);
        double A3D_PowerPattern = std::min((-1) * (AverticalPowerPattern + AhorizontalPowerPattern), 30.0);

        double f1_Theta = sqrt(A3D_PowerPattern) * cos(ksi);
        double f1_Phi = sqrt(A3D_PowerPattern) * sin(ksi);

        fieldPattern << f1_Theta, f1_Phi;

        return fieldPattern;
    }

    //Переход от ЛСК в ГСК 
    Vector2d transformationFromLCSToGCS(double thetaAngle, double phiAngle, double downtiltAngle, Vector2d& fieldPattern) const {
        Vector2d transformFieldPattern;
        double cos_Pci = (cos(downtiltAngle) * sin(thetaAngle) - sin(downtiltAngle) * cos(phiAngle) * cos(thetaAngle)) / (pow((1 - (cos(downtiltAngle) * cos(thetaAngle) - sin(downtiltAngle) * cos(phiAngle) * sin(thetaAngle)) * (cos(downtiltAngle) * cos(thetaAngle) - sin(downtiltAngle) * cos(phiAngle) * sin(thetaAngle))), 0.5));
        double sin_Pci = (sin(downtiltAngle) * sin(phiAngle)) / (pow((1 - (cos(downtiltAngle) * cos(thetaAngle) - sin(downtiltAngle) * cos(phiAngle) * sin(thetaAngle)) * (cos(downtiltAngle) * cos(thetaAngle) - sin(downtiltAngle) * cos(phiAngle) * sin(thetaAngle))), 0.5));

        double f_Theta = sin_Pci * fieldPattern(0) + cos_Pci * fieldPattern(1);
        double f_Phi = cos_Pci * fieldPattern(0) - sin_Pci * fieldPattern(1);

        transformFieldPattern << f_Theta, f_Phi;
        return transformFieldPattern;
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

//_________________________Функция_для_вычисления_расстояния_между_двумя_пользователями____________________//
double calculateDistance(const UserTerminal& a, const UserTerminal& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

//_______________________________________Основная функция программы_________________________________//
int main()
{
    // Параметры комнаты 
    double roomLength = 120.0; // Длина комнаты
    double roomWidth = 50.0;    // Ширина комнаты 
    double roomHeight = 3.0;    // Высота комнаты 
    double wavelength = 0.3;     // Длина волны 
    double nu = 3; // Частота сигнала
    // Создаем пользователей 
    std::vector<UserTerminal> users;
    std::random_device rd;
    std::mt19937_64 gen64(rd());
    std::uniform_real_distribution<> xDist(0.0, roomLength);
    std::uniform_real_distribution<> yDist(0.0, roomWidth);
    double userHeight = 1.0; // Высота пользователей 

    for (int i = 1; i <= 12; ++i) {
        //double bearing = (rand() % 360) * M_PI / 180.0; // Угол поворота 
        double bearing = 0;
        double downtilt = (rand() % 90) * M_PI / 180.0; // Угол наклона
        double slant = 0;
        //double slant = (rand() % 360) * M_PI / 180.0; // Угол наклона

        UserTerminal newUT(i, 0, 0, userHeight, wavelength, bearing, downtilt, slant);
        bool isValidPosition = false;

        while (!isValidPosition) {
            newUT.x = xDist(gen64);
            newUT.y = yDist(gen64);

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
    while (selectedPairs.size() < 6)
    {
        int user1 = userDist(gen64);
        int user2 = userDist(gen64);

        if (user1 != user2) {
            selectedPairs.emplace(std::min(user1, user2), std::max(user1, user2));
        }
    }

    //Генерация PATHLOSS
    // Файлы для Pathloss
    //std::ofstream d_LOS, PATHLOSS_LOS;
    //std::ofstream d_NLOS, PATHLOSS_NLOS;

    //d_LOS.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\LOS\\d_LOS.txt");
    //PATHLOSS_LOS.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\LOS\\PATHLOSS_LOS.txt");

    //d_NLOS.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\NLOS\\d_NLOS.txt");
    //PATHLOSS_NLOS.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\PATHLOSS\\NLOS\\PATHLOSS_NLOS.txt");

    //for (int i = 0; i < 1000; i++)
    //{
    //    // генерация пользователей
    //    for (int j = 1; j <= 12; j++)
    //    {
    //        //double bearing = (rand() % 360) * M_PI / 180.0; // Угол поворота 
    //        double bearing = 0;
    //        double downtilt = (rand() % 90) * M_PI / 180.0; // Угол наклона
    //        double slant = 0;
    //        //double slant = (rand() % 360) * M_PI / 180.0; // Угол наклона

    //        UserTerminal newUT(j, 0, 0, userHeight, wavelength, bearing, downtilt, slant);
    //        bool isValidPosition = false;

    //        while (!isValidPosition) {
    //            newUT.x = xDist(gen64);
    //            newUT.y = yDist(gen64);

    //            // Проверяем расстояние до всех существующих пользователей
    //            isValidPosition = true; // Предполагаем, что позиция валидна
    //            for (const auto& user : users) {
    //                if (calculateDistance(newUT, user) < 1.0) {
    //                    isValidPosition = false; // Позиция невалидна
    //                    break; // Выходим из цикла
    //                }
    //            }
    //            users.push_back(newUT);
    //        }
    //    }
    //    // выбор пар пользователей
    //    std::uniform_int_distribution<> userDist(0, users.size() - 1);
    //    std::set<std::pair<int, int>> selectedPairs;
    //    while (selectedPairs.size() < 6)
    //    {
    //        int user1 = userDist(gen64);
    //        int user2 = userDist(gen64);

    //        if (user1 != user2) {
    //            selectedPairs.emplace(std::min(user1, user2), std::max(user1, user2));
    //        }
    //    }

    //    for (const auto& pair : selectedPairs)
    //    {
    //        const UserTerminal& transmitter = users[pair.first];
    //        const UserTerminal& receiver = users[pair.second];

    //        std::cout << "// UT " << transmitter.id << " - UT " << receiver.id << " // " << std::endl << std::endl;

    //        // Генерация PATHLOSS
    //        //_____________________________________STEP_2________________________________________//

    //        // Вычисление расстояния между передатчиком и приемником
    //        double d = calculateDistance(transmitter, receiver);
    //        std::cout << "Distance between transmitter and receiver is " << d << std::endl;

    //        double P_LOS;
    //        bool los;
    //        if (d <= 5)
    //        {
    //            P_LOS = 1;
    //            los = true;
    //            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
    //        }
    //        else if (d > 5 && d <= 49)
    //        {
    //            P_LOS = exp(-(d - 5) / 70.8);
    //            los = true;
    //            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
    //        }
    //        else
    //        {
    //            P_LOS = exp(-(d - 49) / 211.7) * 0.54;
    //            if (P_LOS < 0.5) {
    //                los = false;
    //                std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - NLOS" << std::endl;
    //            }
    //            else {
    //                los = true;
    //                std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
    //            }
    //        }

    //        //________________________________________STEP_3_____________________________________//
    //        // Вычисление Pathloss_LOS для каждой пары
    //        double nu = 3; // частота в ГГц
    //        if (los)
    //        {
    //            double path_loss_LOS;
    //            path_loss_LOS = 32.4 + 17.3 * log10(d) + 20 * log10(nu);
    //            std::cout << "PathLoss in LOS between users " << transmitter.id << " and " << receiver.id << " = " << path_loss_LOS << std::endl << std::endl;
    //            d_LOS << d << std::endl;
    //            PATHLOSS_LOS << path_loss_LOS << std::endl;
    //        }
    //        else
    //        {
    //            double path_loss_NLOS;
    //            path_loss_NLOS = 38.3 * log10(d) + 17.30 + 24.9 * log10(nu);
    //            std::cout << "Pathloss in NLOS between users " << transmitter.id << " and " << receiver.id << " = " << path_loss_NLOS << std::endl << std::endl;
    //            d_NLOS << d << std::endl;
    //            PATHLOSS_NLOS << path_loss_NLOS << std::endl;
    //        }
    //    }
    //}

    //d_LOS.close();
    //PATHLOSS_LOS.close();

    //d_NLOS.close();
    //PATHLOSS_NLOS.close();


    //// генерация LSP
    for (const auto& pair : selectedPairs)
    {
        const UserTerminal& transmitter = users[pair.first];
        const UserTerminal& receiver = users[pair.second];
        std::cout << "// UT " << transmitter.id << " - UT " << receiver.id << " // " << std::endl << std::endl;
        //________________________________STEP_2___________________________________________________________//
        // Вычисление расстояния между передатчиком и приемником
        double d = calculateDistance(transmitter, receiver);
        std::cout << "Distance between transmitter and receiver is " << d << std::endl;
        double P_LOS;
        bool los;
        if (d <= 5) {
            P_LOS = 1;
            los = true;
            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
        }
        else if (d > 5 && d <= 49) {
            P_LOS = exp(-(d - 5) / 70.8);
            los = true;
            std::cout << "LOS probability between " << transmitter.id << " and " << receiver.id << " = " << P_LOS << ", Link - LOS" << std::endl;
        }
        else {
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
        //_________________________________________________STEP_4__________________________________________//
        std::ofstream SF_LOS_CDF, K_LOS_CDF, DS_LOS_CDF, ASA_LOS_CDF, ASD_LOS_CDF, ZSA_LOS_CDF, ZSD_LOS_CDF;
        std::ofstream SF_LOS_PDF, K_LOS_PDF, DS_LOS_PDF, ASA_LOS_PDF, ASD_LOS_PDF, ZSA_LOS_PDF, ZSD_LOS_PDF;
        std::ofstream SF_NLOS_CDF, DS_NLOS_CDF, ASA_NLOS_CDF, ASD_NLOS_CDF, ZSA_NLOS_CDF, ZSD_NLOS_CDF;
        std::ofstream SF_NLOS_PDF, DS_NLOS_PDF, ASA_NLOS_PDF, ASD_NLOS_PDF, ZSA_NLOS_PDF, ZSD_NLOS_PDF;

    //    // \ - Управляющий символ

        if (los)
        {
    //        // Файлы для CDF
            SF_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\SF_LOS_CDF.txt");
            K_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\K_LOS_CDF.txt");
            DS_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\DS_LOS_CDF.txt");
            ASA_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\ASA_LOS_CDF.txt");
            ASD_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\ASD_LOS_CDF.txt");
            ZSA_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\ZSA_LOS_CDF.txt");
            ZSD_LOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\CDF\\ZSD_LOS_CDF.txt");

    //        // Файлы для PDF
            SF_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\SF_LOS_PDF.txt");
            K_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\K_LOS_PDF.txt");
            DS_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\DS_LOS_PDF.txt");
            ASA_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\ASA_LOS_PDF.txt");
            ASD_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\ASD_LOS_PDF.txt");
            ZSA_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\ZSA_LOS_PDF.txt");
            ZSD_LOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_LOS\\PDF\\ZSD_LOS_PDF.txt");
       }

        else
        {
            // Файлы для CDF
            SF_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\SF_NLOS_CDF.txt");
            DS_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\DS_NLOS_CDF.txt");
            ASA_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\ASA_NLOS_CDF.txt");
            ASD_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\ASD_NLOS_CDF.txt");
            ZSA_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\ZSA_NLOS_CDF.txt");
            ZSD_NLOS_CDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\CDF\\ZSD_NLOS_CDF.txt");

    //        // Файлы для PDF
            SF_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\SF_NLOS_PDF.txt");
            DS_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\DS_NLOS_PDF.txt");
            ASA_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\ASA_NLOS_PDF.txt");
            ASD_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\ASD_NLOS_PDF.txt");
            ZSA_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\ZSA_NLOS_PDF.txt");
            ZSD_NLOS_PDF.open("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 1\\C++\\Graphics\\LSP_NLOS\\PDF\\ZSD_NLOS_PDF.txt");
        }

        for (int i = 0; i < 100; i++)
        {
            // Вывод в зависимости от LOS/NLOS
            if (los)
            {
                std::cout << "\nLSP for LOS for User " << transmitter.id << " and User " << receiver.id << " :\n\n";
            }
            else
            {
                std::cout << "\nLSP for NLOS for User " << transmitter.id << " and User " << receiver.id << " :\n\n";
            }
            LargeScaleParameters lsp(los, nu);
            lsp.showParameters();
            if (los)
            {
             SF_LOS_CDF << lsp.shadowFading << std::endl;
                K_LOS_CDF << lsp.riceanK << std::endl;
                DS_LOS_CDF << lsp.delaySpread << std::endl;
                ASA_LOS_CDF << lsp.azimuthSpreadArrival << std::endl;
                ASD_LOS_CDF << lsp.azimuthSpreadDeparture << std::endl;
                ZSA_LOS_CDF << lsp.zenithSpreadArrival << std::endl;
                ZSD_LOS_CDF << lsp.zenithSpreadDeparture << std::endl;

            SF_LOS_PDF << round(lsp.shadowFading) << std::endl;
                K_LOS_PDF << round(lsp.riceanK * 2) / 2 << std::endl;
                DS_LOS_PDF << round(lsp.delaySpread * 100) / 100 << std::endl;
                ASA_LOS_PDF << round(lsp.azimuthSpreadArrival * 50) / 50 << std::endl;
                ASD_LOS_PDF << round(lsp.azimuthSpreadDeparture * 50) / 50 << std::endl;
                ZSA_LOS_PDF << round(lsp.zenithSpreadArrival * 25) / 25 << std::endl;
                ZSD_LOS_PDF << round(lsp.zenithSpreadDeparture * 50) / 50 << std::endl;
            }

            else
            {
                SF_NLOS_CDF << lsp.shadowFading << std::endl;
                DS_NLOS_CDF << lsp.delaySpread << std::endl;
                ASA_NLOS_CDF << lsp.azimuthSpreadArrival << std::endl;
                ASD_NLOS_CDF << lsp.azimuthSpreadDeparture << std::endl;
                ZSA_NLOS_CDF << lsp.zenithSpreadArrival << std::endl;
                ZSD_NLOS_CDF << lsp.zenithSpreadDeparture << std::endl;

                SF_NLOS_PDF << round(lsp.shadowFading) << std::endl;
                DS_NLOS_PDF << round(lsp.delaySpread * 50) / 50 << std::endl;
                ASA_NLOS_PDF << round(lsp.azimuthSpreadArrival * 50) / 50 << std::endl;
                ASD_NLOS_PDF << round(lsp.azimuthSpreadDeparture * 50) / 50 << std::endl;
                ZSA_NLOS_PDF << round(lsp.zenithSpreadArrival * 100) / 100 << std::endl;
                ZSD_NLOS_PDF << round(lsp.zenithSpreadDeparture * 50) / 50 << std::endl;
           }
        }


        SF_LOS_CDF.close();
        K_LOS_CDF.close();
        DS_LOS_CDF.close();
        ASA_LOS_CDF.close();
        ASD_LOS_CDF.close();
        ZSA_LOS_CDF.close();
        ZSD_LOS_CDF.close();

        SF_LOS_PDF.close();
        K_LOS_PDF.close();
        DS_LOS_PDF.close();
        ASA_LOS_PDF.close();
        ASD_LOS_PDF.close();
        ZSA_LOS_PDF.close();
        ZSD_LOS_PDF.close();

        SF_NLOS_CDF.close();
        DS_NLOS_CDF.close();
        ASA_NLOS_CDF.close();
        ASD_NLOS_CDF.close();
        ZSA_NLOS_CDF.close();
        ZSD_NLOS_CDF.close();

        SF_NLOS_PDF.close();
        DS_NLOS_PDF.close();
        ASA_NLOS_PDF.close();
        ASD_NLOS_PDF.close();
        ZSA_NLOS_PDF.close();
        ZSD_NLOS_PDF.close();

        break;
    }
    return 0;
}