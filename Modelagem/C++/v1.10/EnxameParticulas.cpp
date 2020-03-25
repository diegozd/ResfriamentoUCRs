
#include "EnxameParticulas.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "RESFRIAMENTO.h"
#include "Random.h"

using namespace std;

EnxameParticulas::EnxameParticulas(DADOS &D)
{
    Random R;

    c1 = 2.0;
    c2 = 2.0;
    wIni = 0.9;
    wFim = 0.4;
    vPg = 0.0;
    w = wIni;

    Pg.resize(D.Parametros);
    xMax.resize(D.Parametros);
    xMin.resize(D.Parametros);
    vMax.resize(D.Parametros);

    for (int j = 0 ; j < D.Parametros ; j++)
    {
        xMax[j] = D.xMax[j];
        xMin[j] = D.xMin[j];
        vMax[j] = 0.25*(xMax[j] - xMin[j]);
    }

    FobjValor.resize(D.Particulas);
    vPl.resize(D.Particulas);
    X.resize(D.Particulas);
    V.resize(D.Particulas);
    Pl.resize(D.Particulas);

    for (int i = 0; i < D.Particulas ; i++)
    {
        vPl[i] = 0;
        X[i].resize(D.Parametros);
        V[i].resize(D.Parametros);
        Pl[i].resize(D.Parametros);

        for (int j = 0 ; j < D.Parametros ; j++)
        {
            X[i][j] = R.NumeroRandomico(xMin[j],xMax[j]);
            V[i][j] = vMax[j]*(2*R.NumeroRandomico(0,1)-1);
        }
    }
}


void EnxameParticulas::Otimizar(DADOS &D)
{
    ofstream fout2("PSOout.txt");
	fout2 <<  "Iteracao FuncaoObjetivo" << endl;
	for (int i = 0 ; i < D.Iteracoes; i++)
    {
        AvaliarFobj(D);
        AvaliarSolucao(D);
        w = wFim + ((D.Iteracoes-i)/D.Iteracoes)*(wIni - wFim);
        AtualizarPosicao(D);
		
		fout2 << i << " " << vPg << endl;
		cout <<  "Iteracao " << i << " Funcao Objetivo " << vPg << endl << endl;
    }

    fout2 << endl << "Função objetivo mínima: " << vPg << endl << endl;
    fout2 << "Minimo encontrado em:" << endl;
    for (int i = 0 ; i < D.Parametros ; i++)
    {
        fout2 << "X" << i << " " << Pg[i] << endl;
    }
    fout2 << endl;
    fout2.close();

}

void EnxameParticulas::AvaliarFobj(DADOS &D)
{
    #pragma omp parallel for
	for (int i = 0 ; i < D.Particulas ; i++)
	{
        FobjValor[i] = FuncaoObjetivo(X[i]);
        cout << "Particula: " << i << "  Fobj =  "<< aux << endl;
	}
}

void EnxameParticulas::AvaliarSolucao(DADOS &D)
{
    vector<double>::iterator it;
    vPgN = *min_element(FobjValor.begin(),FobjValor.end());
    it = find(FobjValor.begin(),FobjValor.end(),vPgN);
    int Ind = distance(FobjValor.begin(),it);

    if(vPgN < vPg)
    {
    	//CountConv = 0;
        vPg = vPgN;
        for (int i = 0 ; i < D.Parametros ; i++)
        {
            Pg[i] = X[Ind][i];
        }
    }
    else if (vPg == 0)
    {
    	//CountConv = 0;
        vPg = vPgN;
        for (int i = 0 ; i < D.Parametros ; i++)
        {
            Pg[i] = X[Ind][i];
        }
    }
    else
    {
    	//CountConv += 1;
    }


    for (int i = 0 ; i < D.Particulas ; i++)
    {
        if ( FobjValor[i] < vPl[i] )
        {
            vPl[i] = FobjValor[i];

             for (int j = 0 ; j < D.Parametros ; j++)
             {
                Pl[i][j] = X[i][j];
             }
        }
         else if (vPl[i] == 0)
        {
             vPl[i] = FobjValor[i];
             for (int j = 0 ; j < D.Parametros ; j++)
             {
                Pl[i][j] = X[i][j];
             }
        }
    }
}


void EnxameParticulas::AtualizarPosicao(DADOS &D)
{
    Random R;

	for (int i = 0 ; i < D.Particulas ; i++)
	{
		for(int j = 0 ; j < D.Parametros ; j++)
		{
			V[i][j] = V[i][j]*w + c1*R.NumeroRandomico(0,1)*(Pl[i][j] - X[i][j]) + c2*R.NumeroRandomico(0,1)*(Pg[j] - X[i][j]);

		if (fabs(V[i][j]) > vMax[j])
		{
			if (V[i][j] < 0 )
			{
				V[i][j] = -vMax[j];
			}
			else
			{
				V[i][j] = vMax[j];
			}
		}

			X[i][j] = X[i][j] + V[i][j];


			if (X[i][j] > xMax[j] )
			{
				X[i][j] = xMax[j];
			}
			else if (X[i][j] < xMin[j])
			{
				X[i][j] = xMin[j];
			}
		}
	}
}


double EnxameParticulas::FuncaoObjetivo(vector<double> X)
{
    vector<int> auxPar2;	
	auxPar2.resize(7);
	for (int i = 0 ; i < 7 ; i++)
	{
		auxPar2[i] = 0;
	}
	
	DADOS geral;
	for (int i = 0 ; i < geral.Parametros ; i++)
	{
		if ((geral.auxPar[0] == 1) && (auxPar2[0] == 0))
		{
			geral.Dp = X[i];
			auxPar2[0] = 1;
		}
		else if ((geral.auxPar[1] == 1) && (auxPar2[1] == 0))
		{
			geral.Le = X[i];
			auxPar2[1] = 1;
		}
		else if ((geral.auxPar[2] == 1) && (auxPar2[2] == 0))
		{
			geral.zd = X[i];
			auxPar2[2] = 1;
		}
		else if ((geral.auxPar[3] == 1) && (auxPar2[3] == 0))
		{
			geral.eps = X[i];
			auxPar2[3] = 1;
		}
		else if ((geral.auxPar[4] == 1) && (auxPar2[4] == 0))
		{
			geral.mhc = X[i];
			auxPar2[4] = 1;
		}
		else if((geral.auxPar[5] == 1) && (auxPar2[5] == 0))
		{
			geral.mhl = X[i];
			auxPar2[5] = 1;
		}
		else if((geral.auxPar[6] == 1) && (auxPar2[6] == 0))
		{
			geral.mhe = X[i];
			auxPar2[6] = 1;
		}
	}

    PROPRIEDADES geralPP(geral);
    RESFRIAMENTO geralRESF(geral,geralPP);

    NumDados = geral.ResfExp.size();
    solda = geral.hcone / geral.dL - 1;
    MatrizCompara.resize(NumDados);
    aux = 0;
    int j = 0;
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
            }
            j += 1;
        }while (MatrizCompara[i][0] == -1);
        aux += pow((geral.ResfExp[i][1] - MatrizCompara[i][1]),2)/geral.NddsCrvExp;
    }
	//calculando o erro quadrado para o vapor total gerado no blowdown
	aux += pow((geralRESF.VaporGeradoTotal/1000 - geral.AgVapExp),2);
    return aux;
}
