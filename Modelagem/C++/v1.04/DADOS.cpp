//#include <stdlib.h>
#include "DADOS.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

DADOS::DADOS()
{
	
	//Utilidade p�blica
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
            { 	//Semente da simula��o
                infile >> seed;
            }
            else if (i == 10)
            {	//Di�metro do leito [m]
                infile >> Dc;
            }
			else if (i == 12)
            {	//Altura do Leito de coque [m]
                infile >> AltLeito;
            }
			else if (i == 14)
            {	//Dim�metro menor do cone (di�metro da delta valve) [m]
                infile >> Dcone;
            }
			else if (i == 16)
            {	//Altura do cone (dist�ncia vertical da delta at� a solda)[m]
                infile >> hcone;
            }
			else if (i == 20)
            {	//Press�o do leito durante resfriamento[kgf/cm� man]
                infile >> P;
				P = P / 1.0197 + 1.01325; // [kgf/cm� man >> bar abs] 
            }
			else if (i == 22)
            {	//Temperatura inicial do coque [�C]
                infile >> ToCoq;
				ToCoq = ToCoq + 273.15; // [�C >> K]
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
            {	//Massa espec�fica do coque [kg/m�]
                infile >> rohcoq;
            }
			else if (i == 28+5+k)
            {	//Capacidade calor�fica do coque [kJ/kg.K]
                infile >> cpcoq;
            }
			else if (i == 28+7+k)
            {	//Condutividade t�rmica do coque [W/m.K]
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
            {	//Di�metro m�dio das particulas de coque no leito [m]
                infile >> Dp;
            }
			else if (i == 28+17+k)
            {	//Comprimento da interface/espuma [m]
                infile >> Le;
            }
			else if (i == 28+19+k)
            {	//Par�metro da curva da interface/espuma
                infile >> zd;
            }
			else if (i == 28+21+k)
            {	//Tortuosidade do leito
                infile >> psi;
            }
			else if (i == 28+23+k)
            {	//Par�metro do coeficiente de convec��o for�ada
                infile >> mhc;
            }
			else if (i == 28+25+k)
            {	//Par�metro do coeficiente de convec��o livre
                infile >> mhl;
            }
			else if (i == 28+27+k)
            {	//Par�metro do coeficiente de convec��o em ebuli��o
                infile >> mhe;
            }
			else if (i == 28+31+k)
            {	//Altura do elemento finito de volume [m]
                infile >> dL;
            }
			else if (i == 28+33+k)
            {	//Passo de integra��o no tempo [s]
                infile >> dt;
            }
			else if (i == 28+35+k)
            {	//Numero de amostras de tempo para exibi��o do resultado
                infile >> NamsT;
            }
			else if (i == 28+37+k)
            {	//Numero m�ximo de itera��es por loop
                infile >> nmax;
            }
			else if (i == 28+41+k)
            {	//Se deseja que o seja feito o PSO
                infile >> PSO;
            }
			else if (i == 28+43+k)
            {	//Numero de itera��es do PSO
                infile >> Iteracoes;
            }
			else if (i == 28+45+k)
            {	//Numero de particulas para o PSO
                infile >> Particulas;
            }
			else if (i == 28+47+k)
            {	//Numero de par�metros para o PSO
				Parametros = 0;
				auxPar.resize(7);
				for (int j = 0 ; j < 7 ; j++)
				{
					infile >> auxPar[j];
					Parametros += auxPar[j];
				}
            }
			else if (i == 28+49+k)
            {	//Valores de m�ximo para os par�metros
				xMax.resize(Parametros);
				for (int j = 0 ; j < Parametros ; j++)
				{
					infile >> xMax[j];
				}
            }
			else if (i == 28+51+k)
            {	//Valores de m�nimo para os par�metros
				xMin.resize(Parametros);
				for (int j = 0 ; j < Parametros ; j++)
				{
					infile >> xMin[j];
				}
            }
			else if (i == 28+54+k)
            {	//Total de �gua vaporizada experimental (medido no blowdown)
                infile >> AgVapExp;
            }
			else if (i == 28+56+k)
            {	//Numero de Dados experimentais
                infile >> NddsCrvExp;
				ResfExp.resize(NddsCrvExp);
            }
			else if ((i >= 28+59+k) && (i < 28+59+k+NddsCrvExp))
            {	//Dados experimentais de resfriamento na regi��o de solda
				ResfExp[i-(28+59+k)].resize(2);
				infile >> ResfExp[i-(28+59+k)][0];
				infile >> ResfExp[i-(28+59+k)][1];
            }
		}
		infile.close();
	}
	else
    {
        cout << "arquivo n�o abriu" << endl;
    }
};
	



