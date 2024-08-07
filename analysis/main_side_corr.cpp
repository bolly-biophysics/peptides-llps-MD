#include<iostream>
#include<fstream>
#include<cmath>
#define N_main 40
#define N_side 80
using namespace std;

float sc_center_x[60000][50], sc_center_y[60000][50], sc_center_z[60000][50];
float main_v_x[60000][50], main_v_y[60000][50], main_v_z[60000][50];
float side_v_x[60000][50], side_v_y[60000][50], side_v_z[60000][50];
float mm_corr[50][50][60000], mm_corr_avg[50][50], mm_corr_std[50][50];
float ss_corr[50][50][60000], ss_corr_avg[50][50], ss_corr_std[50][50];
float ms_corr[50][50][60000], ms_corr_avg[50][50], ms_corr_std[50][50];

struct example
{
    char name[5];
    float x, y, z;
};

example main_info[60000][100];
example side_info[60000][100];

int main()
{
    ifstream infile;
    infile.open("main_chain.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将侧链原子坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_main + 0.0); 
        k = fmod(i, N_main);
        if (k == 0) 
        {
            k = N_main;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> main_info[j][k].name >> main_info[j][k].x >> main_info[j][k].y >> main_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    infile.open("side_chain.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将侧链原子坐标格式化存入结构数组
    i = 1;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_side + 0.0); 
        k = fmod(i, N_side);
        if (k == 0) 
        {
            k = N_side;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> side_info[j][k].name >> side_info[j][k].x >> side_info[j][k].y >> side_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile_1, outfile_2, outfile_3, outfile_4;
    outfile_1.open("mm_corr.dat", ios::out | ios::ate);
    outfile_2.open("ss_corr.dat", ios::out | ios::ate);
    outfile_3.open("ms_corr.dat", ios::out | ios::ate);
    outfile_4.open("ms_corr_s.dat", ios::out | ios::ate);

    // 计算侧链电荷中心随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 80; j = j + 2)
        {
            sc_center_x[i][(j + 1) / 2] = (side_info[i][j].x + side_info[i][j + 1].x) / 2;
            sc_center_y[i][(j + 1) / 2] = (side_info[i][j].y + side_info[i][j + 1].y) / 2;
            sc_center_z[i][(j + 1) / 2] = (side_info[i][j].z + side_info[i][j + 1].z) / 2;
        }
    }

    // 计算残基主链部分运动速度矢量的分量
    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            main_v_x[i][j] = main_info[i + 1][j].x - main_info[i][j].x;
            main_v_y[i][j] = main_info[i + 1][j].y - main_info[i][j].y;
            main_v_z[i][j] = main_info[i + 1][j].z - main_info[i][j].z;
        }
    }

    // 计算残基侧链部分运动速度矢量的分量
    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            side_v_x[i][j] = sc_center_x[i + 1][j] - sc_center_x[i][j];
            side_v_y[i][j] = sc_center_y[i + 1][j] - sc_center_y[i][j];
            side_v_z[i][j] = sc_center_z[i + 1][j] - sc_center_z[i][j];
        }
    }
    
    // 计算残基主链-残基主链之间运动速度矢量夹角余弦随时间的变化
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_main; j++)
        {
            for (k = 1; k <= t - 1; k++)
            {
                mm_corr[i][j][k] = (main_v_x[k][i] * main_v_x[k][j] + main_v_y[k][i] * main_v_y[k][j] + main_v_z[k][i] * main_v_z[k][j]) / (sqrt(main_v_x[k][i] * main_v_x[k][i] + main_v_y[k][i] * main_v_y[k][i] + main_v_z[k][i] * main_v_z[k][i]) * sqrt(main_v_x[k][j] * main_v_x[k][j] + main_v_y[k][j] * main_v_y[k][j] + main_v_z[k][j] * main_v_z[k][j]));
            }
        }
    }
    
    // 计算残基主链-残基主链之间运动速度矢量夹角余弦的平均值(方向关联性)
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_main; j++)
        {
            mm_corr_avg[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                mm_corr_avg[i][j] += mm_corr[i][j][k];
            }
            mm_corr_avg[i][j] /= t - 1;
        }
    }
    
    // 计算残基主链-残基主链之间运动速度矢量夹角余弦的标准差
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_main; j++)
        {
            mm_corr_std[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                mm_corr_std[i][j] += (mm_corr[i][j][k] - mm_corr_avg[i][j]) * (mm_corr[i][j][k] - mm_corr_avg[i][j]);
            }
            mm_corr_std[i][j] = sqrt(mm_corr_std[i][j] / (t - 1));
            outfile_1 << i << " " << j << " " << mm_corr_avg[i][j] << " " << mm_corr_std[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();
    
    // 计算残基侧链-残基侧链之间运动速度矢量夹角余弦随时间的变化
    for (i = 1; i <= N_side / 2; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            for (k = 1; k <= t - 1; k++)
            {
                ss_corr[i][j][k] = (side_v_x[k][i] * side_v_x[k][j] + side_v_y[k][i] * side_v_y[k][j] + side_v_z[k][i] * side_v_z[k][j]) / (sqrt(side_v_x[k][i] * side_v_x[k][i] + side_v_y[k][i] * side_v_y[k][i] + side_v_z[k][i] * side_v_z[k][i]) * sqrt(side_v_x[k][j] * side_v_x[k][j] + side_v_y[k][j] * side_v_y[k][j] + side_v_z[k][j] * side_v_z[k][j]));
            }
        }
    }
    
    // 计算残基侧链-残基侧链之间运动速度矢量夹角余弦的平均值(方向关联性)
    for (i = 1; i <= N_side / 2; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            ss_corr_avg[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                ss_corr_avg[i][j] += ss_corr[i][j][k];
            }
            ss_corr_avg[i][j] /= t - 1;
        }
    }
    
    // 计算残基侧链-残基侧链之间运动速度矢量夹角余弦的标准差
    for (i = 1; i <= N_side / 2; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            ss_corr_std[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                ss_corr_std[i][j] += (ss_corr[i][j][k] - ss_corr_avg[i][j]) * (ss_corr[i][j][k] - ss_corr_avg[i][j]);
            }
            ss_corr_std[i][j] = sqrt(ss_corr_std[i][j] / (t - 1));
            outfile_2 << i << " " << j << " " << ss_corr_avg[i][j] << " " << ss_corr_std[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

    // 计算残基主链-残基侧链之间运动速度矢量夹角余弦随时间的变化
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            for (k = 1; k <= t - 1; k++)
            {
                ms_corr[i][j][k] = (main_v_x[k][i] * side_v_x[k][j] + main_v_y[k][i] * side_v_y[k][j] + main_v_z[k][i] * side_v_z[k][j]) / (sqrt(main_v_x[k][i] * main_v_x[k][i] + main_v_y[k][i] * main_v_y[k][i] + main_v_z[k][i] * main_v_z[k][i]) * sqrt(side_v_x[k][j] * side_v_x[k][j] + side_v_y[k][j] * side_v_y[k][j] + side_v_z[k][j] * side_v_z[k][j]));
            }
        }
    }
    
    // 计算残基主链-残基侧链之间运动速度矢量夹角余弦的平均值(方向关联性)
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            ms_corr_avg[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                ms_corr_avg[i][j] += ms_corr[i][j][k];
            }
            ms_corr_avg[i][j] /= t - 1;
        }
    }
    
    // 计算残基主链-残基侧链之间运动速度矢量夹角余弦的标准差
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            ms_corr_std[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                ms_corr_std[i][j] += (ms_corr[i][j][k] - ms_corr_avg[i][j]) * (ms_corr[i][j][k] - ms_corr_avg[i][j]);
            }
            ms_corr_std[i][j] = sqrt(ms_corr_std[i][j] / (t - 1));
            outfile_3 << i << " " << j << " " << ms_corr_avg[i][j] << " " << ms_corr_std[i][j] << endl;
        }
        outfile_3 << endl;
    }
    outfile_3.close();

    // 输出单个残基主链与侧链之间的关联强度
    for (i = 1; i <= N_main; i++)
    {
        for (j = 1; j <= N_side / 2; j++)
        {
            if (i == j)
                outfile_4 << i << " " << ms_corr_avg[i][j] << endl;
        }
    }
    outfile_4.close();

    return 0;
}
