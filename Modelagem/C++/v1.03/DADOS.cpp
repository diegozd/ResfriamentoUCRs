//#include <stdlib.h>
#include "DADOS.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

DADOS::DADOS()
{
	
	//Utilidade pública
    pi = 3.14159265358979323846;

	int k = 0;
	NddsCrvExp = 0;
    ifstream infile("Entrada.txt");
    string strAux;
    if (infile.is_open())
    {
		for (int i = 2 ; i < 28+59+k+NddsCrvExp ; i++)
        {
            getline(infile,strAux);
			if (i == 6)
            { 	//Semente da simulação
                infile >> seed;
            }
            else if (i == 10)
            {	//Diâmetro do leito [m]
                infile >> Dc;
            }
			else if (i == 12)
            {	//Altura do Leito de coque [m]
                infile >> AltLeito;
            }
			else if (i == 14)
            {	//Dimâmetro menor do cone (diâmetro da delta valve) [m]
                infile >> Dcone;
            }
			else if (i == 16)
            {	//Altura do cone (distância vertical da delta até a solda)[m]
                infile >> hcone;
            }
			else if (i == 20)
            {	//Pressão do leito durante resfriamento[kgf/cm² man]
                infile >> P;
				P = P / 1.0197 + 1.01325; // [kgf/cm² man >> bar abs] 
            }
			else if (i == 22)
            {	//Temperatura inicial do coque [°C]
                infile >> ToCoq;
				ToCoq = ToCoq + 273.15; // [°C >> K]
            }
			else if (i == 25)
            {	//Numero de patamares
                infile >> k;
				CurvResf.resize(k);
            }
			else if ((i >= 28) && (i < 28+k))
            {	//Curva de refriamento	
				CurvResf[i-28].resize(3);
				infile >> CurvResf[i-28][0];
				infile >> CurvResf[i-28][1];
				infile >> CurvResf[i-28][2];
            }
			else if (i == 28+3+k)
            {	//Massa específica do coque [kg/m³]
                infile >> rohcoq;
            }
			else if (i == 28+5+k)
            {	//Capacidade calorífica do coque [kJ/kg.K]
                infile >> cpcoq;
            }
			else if (i == 28+7+k)
            {	//Condutividade térmica do coque [W/m.K]
                infile >> Kcoq;
				Kcoq = Kcoq/1000; //[W/m.K >> kW/m.K]
            }
			else if (i == 28+9+k)
            {	//Porosidade do leito
                infile >> eps;
            }
			else if (i == 28+11+k)
            {	//Emissividade do coque
                infile >> emiss;
            }
			else if (i == 28+15+k)
            {	//Diâmetro médio das particulas de coque no leito [m]
                infile >> Dp;
            }
			else if (i == 28+17+k)
            {	//Comprimento da interface/espuma [m]
                infile >> Le;
            }
			else if (i == 28+19+k)
            {	//Parâmetro da curva da interface/espuma
                infile >> zd;
            }
			else if (i == 28+21+k)
            {	//Tortuosidade do leito
                infile >> psi;
            }
			else if (i == 28+23+k)
            {	//Parâmetro do coeficiente de convecção forçada
                infile >> mhc;
            }
			else if (i == 28+25+k)
            {	//Parâmetro do coeficiente de convecção livre
                infile >> mhl;
            }
			else if (i == 28+27+k)
            {	//Parâmetro do coeficiente de convecção em ebulição
                infile >> mhe;
            }
			else if (i == 28+31+k)
            {	//Altura do elemento finito de volume [m]
                infile >> dL;
            }
			else if (i == 28+33+k)
            {	//Passo de integração no tempo [s]
                infile >> dt;
            }
			else if (i == 28+35+k)
            {	//Numero de amostras de tempo para exibição do resultado
                infile >> NamsT;
            }
			else if (i == 28+37+k)
            {	//Numero máximo de iterações por loop
                infile >> nmax;
            }
			else if (i == 28+41+k)
            {	//Se deseja que o seja feito o PSO
                infile >> PSO;
				//cout << "PSO = " << PSO << endl;
            }
			else if (i == 28+43+k)
            {	//Numero de iterações do PSO
                infile >> Iteracoes;
				//cout << "Iteracoes = " << Iteracoes << endl;
            }
			else if (i == 28+45+k)
            {	//Numero de particulas para o PSO
                infile >> Particulas;
            }
			else if (i == 28+47+k)
            {	//Numero de parâmetros para o PSO
                //infile >> Parametros;
				Parametros = 0;
				auxPar.resize(7);
				for (int j = 0 ; j < 7 ; j++)
				{
					infile >> auxPar[j];
					Parametros += auxPar[j];
				}
            }
			else if (i == 28+49+k)
            {	//Valores de máximo para os parâmetros
				xMax.resize(Parametros);
				for (int j = 0 ; j < Parametros ; j++)
				{
					infile >> xMax[j];
				}
            }
			else if (i == 28+51+k)
            {	//Valores de mínimo para os parâmetros
				xMin.resize(Parametros);
				for (int j = 0 ; j < Parametros ; j++)
				{
					infile >> xMin[j];
				}
            }
			else if (i == 28+54+k)
            {	//Total de água vaporizada experimental (medido no blowdown)
                infile >> AgVapExp;
            }
			else if (i == 28+56+k)
            {	//Numero de Dados experimentais
                infile >> NddsCrvExp;
				ResfExp.resize(NddsCrvExp);
				//cout << "NddsCrvExp = " << NddsCrvExp << endl;
            }
			else if ((i >= 28+59+k) && (i < 28+59+k+NddsCrvExp))
            {	//Dados experimentais de resfriamento na regição de solda
				ResfExp[i-(28+59+k)].resize(2);
				infile >> ResfExp[i-(28+59+k)][0];
				infile >> ResfExp[i-(28+59+k)][1];
				
				//cout << i-(28+59+k) << " " <<ResfExp[i-(28+59+k)][0] << " " << ResfExp[i-(28+59+k)][1] << endl;
            }
		}
		infile.close();
	}
	else
    {
        cout << "arquivo não abriu" << endl;
    }
    //system("PAUSE");
};
	
	
    /*
    //DADOS DE GEOMETRIA
    Dc = 6.4;           //diâmetro do leito [m]
    AltLeito = 18.6;    //Altura do Leito de coque [m]
    Dcone = 2.2;        //dimâmetro menor do cone [m]
    hcone = 4.685;      //altura do cone [m]

    //DADOS DE PROCESSO
    P = 3.0 / 1.0197 + 1.01325; // Pressão do leito [bar abs]
    ToCoq = 400.0 + 273.15;     // Temperatura inicial do coque [K]

    CurvResf = {{  3.0,	 50.0, 200.0},  //curva de refriamento
                { 10.0, 140.0, 200.0},  //o numero de linhas é o número de patamares da curva e sem restrições de máximo
                { 12.5, 230.0,  50.0},  //A primeira coluna são as vazões do fluido [t/h]
                { 27.7, 440.0,  50.0},  //A segunda coluna é o tempo  acumulado [minutos]
                { 62.5, 500.0,  50.0},  //A terceira coluna é a temperatura do fluido [°C]
                {222.4, 560.0,  50.0},
                {250.0, 590.0,  50.0}};

    //PROPRIEDADES DO COQUE
    rohcoq = 1350.0;        //massa específica do coque [kg/m³]
    cpcoq = 1.26;           //capacidade calorífica do coque [kJ/kg.K]
    Kcoq = 0.26 / 1000.0;   //condutividade térmica do coque [kW/m.K]
    eps = 0.37;             //porosidade do leito
    emiss = 0.95;           //emissividade do coque

    //PARÂMETROS DE AJUSTE DO MODELO
    psi = 1.0;      //tortuosidade do leito
    Dp = 3.65442;//3.0       //Diâmetro médio das particulas de coque no leito [m]
    mhl = 1.0;      //parâmetro do coeficiente de convecção livre
    mhc = 1.0;      //parâmetro do coeficiente de convecção forçada
    mhe = 1.0;      //parâmetro do coeficiente de convecção ebulição pelicula
    Le = 7.06482;//6.0       //comprimento da interface/espuma [m]
    zd = 0.710938;//0.9       //parâmetro da curva da interface/espuma

    //PARâMETROS NUMÉRICOS DO MODELO
    dL = 0.2;       //Altura do elemento finito de volume [m]
    dt = 20;        //Passo de integração no tempo [s]
    NamsT = 100;    //numero de amostras de tempo para exibição do resultado
    nmax = 100;     //numero máximo de iterações


    //PARA O ENXAME DE PARTÍCULAS
    PSO = true;                 //Se deseja que o seja feito o PSO atriubua 'true' a esta variável, caso não queira atribua 'false'
    Iteracoes = 5;//20             //numero de iterações do PSO
    Particulas = 30;//50            //numero de particulas para o PSO
    Parametros = 3;             //numero de parâmetros para o PSO
    xMax = {6.4, 10.0, 2.0};    //valores de máximo e mínimo para os parâmetros
    xMin = {0.0,  0.0, 0.0};    //Dp, Le, zd

    ResfExp =  {{  1, 401},     //tempo em minutos na 1° coluna e Temperatura na solda em °C na segunda
                { 29, 397},
                { 52, 396},
                { 71, 392},
                { 91, 390},
                {108, 387},
                {125, 384},
                {142, 380},
                {157, 377},
                {172, 373},
                {187, 362},
                {199, 351},
                {212, 337},
                {222, 324},
                {229, 310},
                {240, 292},
                {250, 277},
                {260, 259},
                {273, 241},
                {282, 224},
                {290, 212},
                {299, 200},
                {307, 185},
                {315, 180},
                {323, 165},
                {337, 154},
                {349, 145},
                {360, 140},
                {372, 131},
                {384, 128},
                {398, 121},
                {410, 117},
                {418, 115},
                {429, 114},
                {444, 109},
                {457, 103},
                {470,  97},
                {481,  93},
                {493,  89},
                {505,  85},
                {517,  83},
                {531,  80},
                {540,  77},
                {550,  75}};
	*/


