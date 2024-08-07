#include<iostream>
#include<fstream>
#include<cmath>
#define bin 33
#define KB 0.001987
#define T 300
using namespace std;

double a[100000], b[100000], c[100000];
double P_original[1000][1000], P_reweight[1000][1000];
double pmf_original[1000][1000], pmf_reweight[1000][1000];
double Boltzmann_factor[1000][1000] = {0};

int main()
{
    ifstream infile;
    infile.open("sb_sasa_weight.dat");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将原子坐标格式化存入结构数组
    int i = 1, j, k, t, count;
    double aa;
    while (!infile.eof())
    {
        infile >> a[i] >> b[i] >> c[i];
        i++;
    }
    t = i - 2;
    cout << t << endl;
    infile.close();

    ofstream outfile_1, outfile_2;
    outfile_1.open("pmf_original.dat", ios::out | ios::ate);
    outfile_2.open("pmf_reweight.dat", ios::out | ios::ate);

    // 计算a列上下限并打印
    double a_max = 0, a_min = 10000;
    for (i = 1; i <= t; i++)
    {
        if (a[i] > a_max)
            a_max = a[i];
        else if (a[i] < a_min)
            a_min = a[i];
    }
    cout << a_min << " " << a_max << endl;

    // 计算b列上下限并打印
    double b_max = 0, b_min = 10000;
    for (i = 1; i <= t; i++)
    {
        if (b[i] > b_max)
            b_max = b[i];
        else if (b[i] < b_min)
            b_min = b[i];
    }
    cout << b_min << " " << b_max << endl;

    // 计算a列、b列的区块跨度
    double a_bin_span = (a_max - a_min) / bin;
    double b_bin_span = (b_max - b_min) / bin;

    // 计算原始二维概率分布
    double pmf_original_min = 100;
    for (i = 0; i <= bin - 1; i++)
    {
        for (j = 0; j <= bin - 1; j++)
        {
            count = 0;
            for (k = 1; k <= t; k++)
            {
                    if (a[k] >= a_min + i * a_bin_span && a[k] < a_min + (i + 1) * a_bin_span)
                    {
                        if (b[k] >= b_min + j * b_bin_span && b[k] < b_min + (j + 1) * b_bin_span)
                            count++;
                    }
            }
            P_original[i][j] = (double)count / t;
            pmf_original[i][j] = -(KB * T) * log(P_original[i][j]);
            if (pmf_original[i][j] < pmf_original_min)
                pmf_original_min = pmf_original[i][j];
        }
    }

    // 输出原始PMF
    for (i = 0; i <= bin - 1; i++)
    {
        for (j = 0; j <= bin - 1; j++)
        {
            outfile_1 << a_min + i * a_bin_span << " " << b_min + j * b_bin_span << " " << pmf_original[i][j] - pmf_original_min << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();

    // 对每个区块计算Ensemble-averaged Boltzmann factor
    double Boltzmann_factor_sum = 0;
    for (i = 0; i <= bin - 1; i++)
    {
        for (j = 0; j <= bin - 1; j++)
        {
            count = 0;
            for (k = 1; k <= t; k++)
            {
                    if (a[k] >= a_min + i * a_bin_span && a[k] < a_min + (i + 1) * a_bin_span)
                    {
                        if (b[k] >= b_min + j * b_bin_span && b[k] < b_min + (j + 1) * b_bin_span)
                        {
                            Boltzmann_factor[i][j] += exp(c[k] / (KB * T));
                            count++;
                        }
                    }
            }
            if (count == 0)
                continue;
            else
                Boltzmann_factor[i][j] = Boltzmann_factor[i][j] / (double)count;
            Boltzmann_factor_sum += P_original[i][j] * Boltzmann_factor[i][j];
        }
    }

    // 计算加权二维概率分布
    double pmf_reweight_min = 100;
    for (i = 0; i <= bin - 1; i++)
    {
        for (j = 0; j <= bin - 1; j++)
        {
            P_reweight[i][j] = P_original[i][j] * Boltzmann_factor[i][j] / Boltzmann_factor_sum;
            pmf_reweight[i][j] = -(KB * T) * log(P_reweight[i][j]);
            if (pmf_reweight[i][j] < pmf_reweight_min)
                pmf_reweight_min = pmf_reweight[i][j];
        }
    }

    // 输出加权PMF
    for (i = 0; i <= bin - 1; i++)
    {
        for (j = 0; j <= bin - 1; j++)
        {
            if ((pmf_reweight[i][j] - pmf_original_min) > 50)
            outfile_2 << a_min + i * a_bin_span << " " << b_min + j * b_bin_span << " " << "inf" << endl;
            else
            outfile_2 << a_min + i * a_bin_span << " " << b_min + j * b_bin_span << " " << pmf_reweight[i][j] - pmf_reweight_min << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();
    
    return 0;
}
