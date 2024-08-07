#include<iostream>
#include<fstream>
#include<cmath>
#define N_ATOM 764
using namespace std;

float d_NA1_O1[10000], d_NA1_O2[10000], d_NA1_O3[10000], d_NA1_O4[10000], d_NA2_O1[10000], d_NA2_O2[10000], d_NA2_O3[10000], d_NA2_O4[10000];

struct example
{
    char name[5];
    float x, y, z;
};

example atom_info[10000][800];

int main()
{
    ifstream infile;
    infile.open("traj.xyz");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将P原子坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_ATOM + 0.0); 
        k = fmod(i, N_ATOM);
        if (k == 0) 
        {
            k = N_ATOM;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> atom_info[j][k].name >> atom_info[j][k].x >> atom_info[j][k].y >> atom_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile;
    outfile.open("dist_O_Na.dat", ios::out | ios::ate);

    // 计算Na离子与O之间的距离
    for (i = 1; i <= t; i++)
    {
        d_NA1_O1[i] = sqrt((atom_info[i][43].x - atom_info[i][13].x) * (atom_info[i][43].x - atom_info[i][13].x) + (atom_info[i][43].y - atom_info[i][13].y) * (atom_info[i][43].y - atom_info[i][13].y) + (atom_info[i][43].z - atom_info[i][13].z) * (atom_info[i][43].z - atom_info[i][13].z));
        d_NA1_O2[i] = sqrt((atom_info[i][43].x - atom_info[i][14].x) * (atom_info[i][43].x - atom_info[i][14].x) + (atom_info[i][43].y - atom_info[i][14].y) * (atom_info[i][43].y - atom_info[i][14].y) + (atom_info[i][43].z - atom_info[i][14].z) * (atom_info[i][43].z - atom_info[i][14].z));
        d_NA1_O3[i] = sqrt((atom_info[i][43].x - atom_info[i][34].x) * (atom_info[i][43].x - atom_info[i][34].x) + (atom_info[i][43].y - atom_info[i][34].y) * (atom_info[i][43].y - atom_info[i][34].y) + (atom_info[i][43].z - atom_info[i][34].z) * (atom_info[i][43].z - atom_info[i][34].z));
        d_NA1_O4[i] = sqrt((atom_info[i][43].x - atom_info[i][35].x) * (atom_info[i][43].x - atom_info[i][35].x) + (atom_info[i][43].y - atom_info[i][35].y) * (atom_info[i][43].y - atom_info[i][35].y) + (atom_info[i][43].z - atom_info[i][35].z) * (atom_info[i][43].z - atom_info[i][35].z));
        d_NA2_O1[i] = sqrt((atom_info[i][44].x - atom_info[i][13].x) * (atom_info[i][44].x - atom_info[i][13].x) + (atom_info[i][44].y - atom_info[i][13].y) * (atom_info[i][44].y - atom_info[i][13].y) + (atom_info[i][44].z - atom_info[i][13].z) * (atom_info[i][44].z - atom_info[i][13].z));
        d_NA2_O2[i] = sqrt((atom_info[i][44].x - atom_info[i][14].x) * (atom_info[i][44].x - atom_info[i][14].x) + (atom_info[i][44].y - atom_info[i][14].y) * (atom_info[i][44].y - atom_info[i][14].y) + (atom_info[i][44].z - atom_info[i][14].z) * (atom_info[i][44].z - atom_info[i][14].z));
        d_NA2_O3[i] = sqrt((atom_info[i][44].x - atom_info[i][34].x) * (atom_info[i][44].x - atom_info[i][34].x) + (atom_info[i][44].y - atom_info[i][34].y) * (atom_info[i][44].y - atom_info[i][34].y) + (atom_info[i][44].z - atom_info[i][34].z) * (atom_info[i][44].z - atom_info[i][34].z));
        d_NA2_O4[i] = sqrt((atom_info[i][44].x - atom_info[i][35].x) * (atom_info[i][44].x - atom_info[i][35].x) + (atom_info[i][44].y - atom_info[i][35].y) * (atom_info[i][44].y - atom_info[i][35].y) + (atom_info[i][44].z - atom_info[i][35].z) * (atom_info[i][44].z - atom_info[i][35].z));
    }

    // 输出八对距离
    for (i = 1; i <= t; i++)
    {
        outfile << i * 2 << " " << d_NA1_O1[i] << " " << d_NA1_O2[i] << " " << d_NA1_O3[i] << " " << d_NA1_O4[i] << " " << d_NA2_O1[i] << " " << d_NA2_O2[i] << " " << d_NA2_O3[i] << " " << d_NA2_O4[i] << endl;
    }
    outfile.close();

    return 0;
}
