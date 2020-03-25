#ifndef CONDUCAO_H
#define CONDUCAO_H

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "CONVECCAO.h"

using namespace std;

class CONDUCAO
{
    public:

        string MsgConducao;

        CONDUCAO();

        double CalcTsup(DADOS&, PROPRIEDADES&, CONVECCAO&);
        double CalcTcoq(DADOS&, PROPRIEDADES&, CONVECCAO&);
        double raizksi(double, int, DADOS&);


    private:

        double Tsuper;      //Tempertaura da superfície dado elemento de volume no proximo intervalo de tempo [K]
        double Tcok;        //Tempertaura do coque num dado elemento de volume no próximo dado intervalo de tempo [K]
        double Tcoq0;       //Tempertaura inicial do coque num dado elemento de volume num dado intervalo de tempo [K]
        double Bi;          //Número de Biot
        double alphaCoq;    //difusividade térmica do coque [m²/s]
        double Fo;          //Número de Fourrier
        double tetaold;     //auxiliar para somatório de teta
        double teta;        //somatório de teta
        double ksi;         //raiz n da função 1-ksi*cot(ksi)=Bi
        double dksi;        //delta de ksi para o cálculo da derivada de fksi
        double ksidksi;     //ksi acrecido do delta de ksi  para o cálculo da derivada de fksi
        double fksi;        //Função Residuo de ksi
        double fksidksi;    //Função Residuo de ksi acrecido do delta de ksi  para o cálculo da derivada de fksi
        double flksi;       //derivada da Função Residuo de ksi
        double aux;         //Variável auxiliar para o cálculo de ksi
        double Cn;          //Variável Cn do modelo de resfriamneto em esfera
        double DIFF1;       //Função Residuo da convergência da serie
        double QQ;          //calor roubado daquele elemento de volume no tempo dt [kJ]
        double DTcoq;       //DTml do coque neste elemento de volume [K]
        double dttt;        //delta de T para o cálculo da derivada de fT
        double DT1;         //variável auxiliar para o cálculo do DTml
        double DT2;         //variável auxiliar para o cálculo do DTml
        double Tsup0;       //Tempertaura inicial na superfície do coque num dado elemento de volume [K]
        double fT;          //Função Residuo de T
        double fTdT;        //Função Residuo de T acrecido do delta de T  para o cálculo da derivada de fT
        double flT;         //derivada da Função Residuo de T

        int n;              //variável que conta o número de termos da série
		//int aviso;
		int aviso2;
		int aviso3;
		int aviso4;
};

#endif // CONDUCAO_H

       // virtual ~CONDUCAO();

   // protected:

   // private:
