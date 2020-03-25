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

        double Tsuper;      //Tempertaura da superf�cie dado elemento de volume no proximo intervalo de tempo [K]
        double Tcok;        //Tempertaura do coque num dado elemento de volume no pr�ximo dado intervalo de tempo [K]
        double Tcoq0;       //Tempertaura inicial do coque num dado elemento de volume num dado intervalo de tempo [K]
        double Bi;          //N�mero de Biot
        double alphaCoq;    //difusividade t�rmica do coque [m�/s]
        double Fo;          //N�mero de Fourrier
        double tetaold;     //auxiliar para somat�rio de teta
        double teta;        //somat�rio de teta
        double ksi;         //raiz n da fun��o 1-ksi*cot(ksi)=Bi
        double dksi;        //delta de ksi para o c�lculo da derivada de fksi
        double ksidksi;     //ksi acrecido do delta de ksi  para o c�lculo da derivada de fksi
        double fksi;        //Fun��o Residuo de ksi
        double fksidksi;    //Fun��o Residuo de ksi acrecido do delta de ksi  para o c�lculo da derivada de fksi
        double flksi;       //derivada da Fun��o Residuo de ksi
        double aux;         //Vari�vel auxiliar para o c�lculo de ksi
        double Cn;          //Vari�vel Cn do modelo de resfriamneto em esfera
        double DIFF1;       //Fun��o Residuo da converg�ncia da serie
        double QQ;          //calor roubado daquele elemento de volume no tempo dt [kJ]
        double DTcoq;       //DTml do coque neste elemento de volume [K]
        double dttt;        //delta de T para o c�lculo da derivada de fT
        double DT1;         //vari�vel auxiliar para o c�lculo do DTml
        double DT2;         //vari�vel auxiliar para o c�lculo do DTml
        double Tsup0;       //Tempertaura inicial na superf�cie do coque num dado elemento de volume [K]
        double fT;          //Fun��o Residuo de T
        double fTdT;        //Fun��o Residuo de T acrecido do delta de T  para o c�lculo da derivada de fT
        double flT;         //derivada da Fun��o Residuo de T

        int n;              //vari�vel que conta o n�mero de termos da s�rie
		//int aviso;
		int aviso2;
		int aviso3;
		int aviso4;
};

#endif // CONDUCAO_H

       // virtual ~CONDUCAO();

   // protected:

   // private:
