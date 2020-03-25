#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "CONVECCAO.h"
#include "CONDUCAO.h"
#include "RESFRIAMENTO.h"
#include "Random.h"
#include "EnxameParticulas.h"

using namespace std;

double FuncaoObjetivoSimples(DADOS &geral, double VaporGeradoTotal)
{
    vector<vector<double> > MatrizCompara;
    int NumDados = geral.NddsCrvExp;
    double solda = geral.hcone / geral.dL - 1;
    double aux = 0;
    int j = 0;
    MatrizCompara.resize(NumDados);
		
    ofstream fout1("CurvaFuncObj.txt");
    fout1 << "t(min)" << "   " << "TckSld(°C)" << endl;
    for (int i = 0 ; i < NumDados ; i++)
    {	//calculando o erro quadrado para os pontos da curva de Temperatura na solda ponderado pelo numero de pontos
        MatrizCompara[i].resize(2);
        MatrizCompara[i][0] = -1;
        do
        {		
            if (j*geral.dt/60 >= geral.ResfExp[i][0])
            {
                MatrizCompara[i][0] =  j*geral.dt/60;
                MatrizCompara[i][1] =  geral.Tcoq[j][solda] - 273.15;
                fout1 << MatrizCompara[i][0] << "   " << MatrizCompara[i][1] << endl;
            }
            j += 1;
        }while (MatrizCompara[i][0] == -1);

        aux += pow((geral.ResfExp[i][1] - MatrizCompara[i][1]),2)/geral.NddsCrvExp;
		
    }
	
	aux += pow((VaporGeradoTotal/1000 - geral.AgVapExp),2);
	
    fout1 << endl << "FuncObj: " << aux << endl;
    fout1.close();
	
    return aux;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

time_t start_time,end_time;

int main()
{
    start_time = time(NULL);
    
	//Informações de abertura
	cout << "PETROBRAS" << endl;
	cout << "Modelo de Resfriamento do Tambor de Coque" << endl;
	cout << "Diego Telles (U3YP)" << endl;
	cout << "Leonardo Caliman (CDW7)" << endl;
	cout << "Versao 1.03 - 05/12/2018" << endl;
    cout << endl << "Calma que ta rodando!" << endl << endl;
	
	DADOS D;
	
	if (D.PSO == 0)
	{
		cout << "Warnings:" << endl;
	}
	else
	{
		cout << "Realizando exame de particulas" << endl;
	}
	

    //realizando calculos de resfirmaneto
    
    srand(D.seed); 
	
	vector<int> auxPar2;
	auxPar2.resize(7);
	for (int i = 0 ; i < 7 ; i++)
	{
		auxPar2[i] = 0;
	}
    if (D.PSO == 1)
    {
        EnxameParticulas EP(D);
        EP.Otimizar(D);
		for (int i = 0 ; i < D.Parametros ; i++)
		{
			if ((D.auxPar[0] == 1) && (auxPar2[0] == 0))
			{
				D.Dp = EP.Pg[i];
				auxPar2[0] = 1;
			}
			else if ((D.auxPar[1] == 1) && (auxPar2[1] == 0))
			{
				D.Le = EP.Pg[i];
				auxPar2[1] = 1;
			}
			else if ((D.auxPar[2] == 1) && (auxPar2[2] == 0))
			{
				D.zd = EP.Pg[i];
				auxPar2[2] = 1;
			}
			else if ((D.auxPar[3] == 1) && (auxPar2[3] == 0))
			{
				D.eps = EP.Pg[i];
				auxPar2[3] = 1;
			}
			else if ((D.auxPar[4] == 1) && (auxPar2[4] == 0))
			{
				D.mhc = EP.Pg[i];
				auxPar2[4] = 1;
			}
			else if((D.auxPar[5] == 1) && (auxPar2[5] == 0))
			{
				D.mhl = EP.Pg[i];
				auxPar2[5] = 1;
			}
			else if((D.auxPar[6] == 1) && (auxPar2[6] == 0))
			{
				D.mhe = EP.Pg[i];
				auxPar2[6] = 1;
			}
		}
    }
	
    PROPRIEDADES PP(D);
    RESFRIAMENTO RESF(D,PP);
	
	double minFobj = 0;
	if (D.PSO == 1) minFobj = FuncaoObjetivoSimples(D, RESF.VaporGeradoTotal);


    //declarando variáveis para os resultados
    int dTT;
    int solda;
    int X;
    int aux4;
    double aux1;
    double aux2;
    double aux3;
    double balanco;
    double errobalanco;
    vector<vector<double> > RESULTADOS;  //Matriz de resultados
	
	solda = D.hcone / D.dL - 1;
	//escrevendo o vetor de derivadas
	ofstream fout3("derivadas.txt");
	//fout3 << "Curva de derivadas da temperatura da altura da solda" << endl;
	fout3 << "tempo(min) dT/dt(°C/min)" << endl;
	for (int i = 1 ; i < D.NumT ; i++)
	{
		aux1 = (D.Tcoq[i-1][solda] - D.Tcoq[i][solda])*60/D.dt;     
		fout3 << i*D.dt/60 << " " << aux1 << endl;
	}
	fout3.close();
	
	//escrevendo resultados finais do programa
    dTT = D.NumT / D.NamsT;
    aux1 = D.NumT;
    aux2 = D.NamsT;
    aux3 = aux1/aux2;
    if ((aux3 - dTT) > 0)
    {
        aux1 = D.NumT - dTT * D.NamsT;
        aux4 = aux1/dTT;
        aux4 += D.NamsT;
    }
    else
    {
        aux4 = D.NamsT;
    }
    RESULTADOS.resize(aux4+1);
    
    X = D.AltLeito/D.dL - 1;
    for(int i = 0 ; i < aux4 ; i++)
    {
        RESULTADOS[i].resize(9);

        RESULTADOS[i][0] = (i * dTT * D.dt) / 60;               //tempo da amostra [min]
        RESULTADOS[i][1] = D.Tflu[i * dTT][D.NumL] - 273.15; 	//Temperatura do vapor na saída do leito [°C]
        RESULTADOS[i][2] = D.Tcoq[i * dTT][0] - 273.15;         //Temperatura do coque no início do leito [°C]
        RESULTADOS[i][3] = D.Tcoq[i * dTT][solda] - 273.15;     //Temperatura do coque na região da solda [°C]
        RESULTADOS[i][4] = D.Tcoq[i * dTT][X] - 273.15;         //Temperatura do coque no topo do leito [°C]
        RESULTADOS[i][5] = D.mvap[i * dTT][D.NumL] * 3.6;    	//vazão de vapor que sai do leito [t/h]
        RESULTADOS[i][6] = RESF.NIVEL[i * dTT];                 //Nível [m]
        RESULTADOS[i][7] = D.Tsup[i * dTT][0] - 273.15;         //Temperatura da superficie do coque no início do leito [°C]
        RESULTADOS[i][8] = D.Tsup[i * dTT][solda] - 273.15;     //Temperatura da superficie do coque na região da solda [°C]

    }

    RESULTADOS[aux4].resize(9);

    RESULTADOS[aux4][0] = D.CurvResf[(D.CurvResf.size()-1)][1];	//tempo da amostra [min]
    RESULTADOS[aux4][1] = D.Tflu[D.NumT-1][D.NumL] - 273.15;	//Temperatura do vapor na saída do leito [°C]
    RESULTADOS[aux4][2] = D.Tcoq[D.NumT][0] - 273.15; 			//Temperatura do coque no início do leito [°C]
    RESULTADOS[aux4][3] = D.Tcoq[D.NumT][solda] - 273.15;   	//Temperatura do coque na região da solda [°C]
    RESULTADOS[aux4][4] = D.Tcoq[D.NumT][X] - 273.15;     		//Temperatura do coque no topo do leito [°C]
    RESULTADOS[aux4][5] = D.mvap[D.NumT-1][D.NumL] * 3.6;  		//vazão de vapor que sai do leito [t/h]
    RESULTADOS[aux4][6] = RESF.NIVEL[D.NumT];              		//Nível [m]
    RESULTADOS[aux4][7] = D.Tsup[D.NumT][0] - 273.15;        	//Temperatura da superficie do coque no início do leito [°C]
    RESULTADOS[aux4][8] = D.Tsup[D.NumT][solda] - 273.15;   	//Temperatura da superficie do coque na região da solda [°C]

    balanco = RESF.TotalAgua - (RESF.VaporGeradoTotal + RESF.AgAcumu);
    errobalanco = 100*balanco/RESF.TotalAgua;

    //Escrevendo resultados
    cout << endl << "Escrevendo Resultados" << endl << endl;
    ofstream fout("Saida.txt");
    fout << "t(min)" << "   " << "Tvap(°C)" << "   " << "TckIn(°C)" << "   " << "TckSld(°C)" << "   " << "TckOut(°C)" << "   " << "Mvap(t/h)" << "   " << "Nivel(m)" << "   " << "TspIn(°C)" << "   " << "TspSld(°C)" << endl;
    for(int i = 0 ; i < aux4+1 ; i++)
    {
        for(int j = 0 ; j < 9 ; j++)
        {
            fout << RESULTADOS[i][j] << "   ";
        }
        fout << endl;
    }
	
    fout << endl;
	if (D.PSO == 1)
    {
		fout << "Função objetivo mínima:"<< " " << minFobj <<endl;
		fout << "Parâmetros otimizados:" << endl;
		for (int i = 0 ; i < 7 ; i++)
		{
			auxPar2[i] = 0;
		}
		for (int i = 0 ; i < D.Parametros ; i++)
		{
			if ((D.auxPar[0] == 1) && (auxPar2[0] == 0))
			{
				fout << "Dp = " << D.Dp << endl;
				auxPar2[0] = 1;
			}
			else if ((D.auxPar[1] == 1) && (auxPar2[1] == 0))
			{
				fout << "Le = " << D.Le << endl;
				auxPar2[1] = 1;
			}
			else if ((D.auxPar[2] == 1) && (auxPar2[2] == 0))
			{
				fout << "zd = " << D.zd << endl;
				auxPar2[2] = 1;
			}
			else if ((D.auxPar[3] == 1) && (auxPar2[3] == 0))
			{
				fout << "eps = " << D.eps << endl;
				auxPar2[3] = 1;
			}
			else if ((D.auxPar[4] == 1) && (auxPar2[4] == 0))
			{
				fout << "mhc = " << D.mhc << endl;
				auxPar2[4] = 1;
			}
			else if((D.auxPar[5] == 1) && (auxPar2[5] == 0))
			{
				fout << "mhl = " << D.mhl << endl;
				auxPar2[5] = 1;
			}
			else if((D.auxPar[6] == 1) && (auxPar2[6] == 0))
			{
				fout << "mhe = " << D.mhe << endl;
				auxPar2[6] = 1;
			}
		}
		fout << "Numero de iterações = " << D.Iteracoes << endl;
		fout << "Numero de particulas = " << D.Particulas << endl;
	}
    fout << endl << "Pico de vazao de vapor [t/h] = " << RESF.PicoVazao*3.6 << endl;	//Pico de vazão de vapor [kg/s]
    fout << "Agua total adicionada [kg] = " << RESF.TotalAgua << endl;            		//totalizador da água que foi injetada no leito [kg]
    fout << "Vapor total gerado [kg] = " << RESF.VaporGeradoTotal << endl;       		//totalizador de todo vapor gerado no leito [kg]
    fout << "Agua acumulada no leito [kg] = " << RESF.AgAcumu << endl;          		//totalizador da água acumulada no leito [kg]
    fout << "balanco [kg] = " << balanco << endl;
    fout << "errobalanco [%] = " << errobalanco << endl << endl;

    end_time = time(NULL);

    fout << "Tempo: "<< end_time-start_time << " s"<< endl;
    fout.close();

	
    int  auxtempo = 0;
	int  auxtempo2 = 0;
	int  auxtempo3 = 0;
	int TempoInteiro;
	
	TempoInteiro = end_time-start_time;
	if (TempoInteiro > 3600)
	{
		auxtempo = TempoInteiro/3600;
		auxtempo3 = TempoInteiro - auxtempo*3600;
		auxtempo2 = auxtempo3/60;
		auxtempo3 = TempoInteiro - auxtempo*3600 - auxtempo2*60;
		cout << "Tempo de execucao: " << auxtempo << "h " << auxtempo2 << "min " << auxtempo3 << "s." << endl;
	}
	else if (TempoInteiro > 60)
	{
		auxtempo = TempoInteiro/60;
		auxtempo2 = TempoInteiro - auxtempo*60;
		cout << "Tempo de execucao: " << auxtempo << "min " << auxtempo2 << "s." << endl;
		
	}
	else
	{
		cout << "Tempo de execucao: " << TempoInteiro << "s." << endl;
	}
	
	cout << endl;
	cout << "FFFFFFFFFFFF   IIII   MMMM     MMMM" << endl;
	cout << "FFFFFFFFFFFF   IIII   MMMMM   MMMMM" << endl;
	cout << "FFFFFFFFFFFF   IIII   MMMMMM MMMMMM" << endl;
	cout << "FFFF                  MMMMMMMMMMMMM" << endl;
	cout << "FFFF           IIII   MMMM MMM MMMM" << endl;
	cout << "FFFFFFFFF      IIII   MMMM  M  MMMM" << endl;
	cout << "FFFFFFFFF      IIII   MMMM     MMMM" << endl;
	cout << "FFF            IIII   MMMM     MMMM" << endl;
	cout << "FFF            IIII   MMMM     MMMM" << endl;
	cout << "FFF            IIII   MMMM     MMMM" << endl;
	cout << endl;
	
	system("PAUSE");

    return 0;
}

