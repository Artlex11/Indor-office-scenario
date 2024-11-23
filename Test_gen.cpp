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

            GaussGenerator rand(0, 1);
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


int main()
{
    std::vector<double> values; // Вектор для хранения сгенерированных значений
    GaussGenerator rand(0, 1);
    
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

    bool los = true;
    if (los)
    {
        //        // Файлы для CDF
        SF_LOS_CDF.open("SF_LOS_CDF.txt");
        K_LOS_CDF.open("K_LOS_CDF.txt");
        DS_LOS_CDF.open("DS_LOS_CDF.txt");
        ASA_LOS_CDF.open("ASA_LOS_CDF.txt");
        ASD_LOS_CDF.open("ASD_LOS_CDF.txt");
        ZSA_LOS_CDF.open("ZSA_LOS_CDF.txt");
        ZSD_LOS_CDF.open("ZSD_LOS_CDF.txt");


    }

    else
    {
        // Файлы для CDF
        SF_NLOS_CDF.open("SF_NLOS_CDF.txt");
        DS_NLOS_CDF.open("DS_NLOS_CDF.txt");
        ASA_NLOS_CDF.open("ASA_NLOS_CDF.txt");
        ASD_NLOS_CDF.open("ASD_NLOS_CDF.txt");
        ZSA_NLOS_CDF.open("ZSA_NLOS_CDF.txt");
        ZSD_NLOS_CDF.open("ZSD_NLOS_CDF.txt");


    }

    for (int i = 0; i < 100000; i++)
    {

        LargeScaleParameters lsp(los, 30);
        values.push_back(lsp.shadowFading); // Добавление сгенерированного значения в вектор

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
    }


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


    std::sort(values.begin(), values.end());

    int unique_count = 0;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i == 0 || values[i] != values[i - 1]) {
            unique_count++;
        }
    }



    // Вывод количества уникальных значений
    std::cout << "unique count: " << unique_count << std::endl;

    std::ifstream file("C:\\Users\\RadioChelik322\\source\\repos\\Test_gen\\SF_LOS_CDF.txt"); // Замените 'yourfile.txt' на имя вашего файла
    if (!file) {
        std::cerr << "Ошибка при открытии файла!" << std::endl;
        return 1; // Завершаем программу, если файл не удалось открыть
    }

    std::set<std::string> unique_values; // Создаем множество для хранения уникальных значений
    std::string value;

    // Читаем значения из файла
    while (file >> value) {
        unique_values.insert(value); // Добавляем значение в множество (дубликаты игнорируются)
    }

    // Выводим количество уникальных значений
    std::cout << "Количество уникальных значений: " << unique_values.size() << std::endl;

    file.close(); // Закрываем файл

    
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