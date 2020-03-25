#ifndef CONVECCAO_H
#define CONVECCAO_H

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"
#include "PROPRIEDADES.h"

using namespace std;

class CONVECCAO
{
    public:

        double Tfilme;  //Temperatura do filme [K]
        double Tbulk;   //Temperatura do bulk [K]
        double h;       //coeficiente de convec��o combindo (livre + for�ada em meio poroso) [kW/m�.K]
        double q;       //Calor de convec��o convergido [kJ/s]
		double Qmaximo; //Calor m�ximo que pode ser trocado pelo fluido.

        string MsgConveccao;

        CONVECCAO();

        double CalcTsai(DADOS&, PROPRIEDADES&, int);
        double calch(double, double, double, double, DADOS&, PROPRIEDADES&, int);
        double calchebu(double, DADOS&, PROPRIEDADES&);
        double calchmax(DADOS&, PROPRIEDADES&);


    private:

        double T0;      //Temperatura de entrada no elemeto finito [K]
        double Tsup;    //Temperatura m�dia da superf�cio num dado elemeto finito [K]
        double Tout;    //Temperatura de sa�da do elemnto para converg�ncia [K]
        double TSai;    //Temperatura de sa�da auxoiliar  do elemento de volume [K]
        double TSail;   //Temepratura de sa�da + incremento de Temperatura para o c�lculo da derivada [K]
        double Tbulkl;  //Temperatura do bulk + incremento de Temperatura para o c�lculo da derivada [K]
        double dTT;     //incremento de temperatura [K]
        double DT1;     //vari�vel auxiliar para o c�lculo do DTml
        double DT2;     //vari�vel auxiliar para o c�lculo do DTml
        double DTml;    //Delta T m�dio logar�timo [K]
        double DTmll;   //Delta T m�dio logar�timo levando em conta o incremento de T (dTT) para o c�lculo da derivada [K]

        double Hin;     //Entalpia do fluido na estrada do elemento de volume [kJ/kg]
        double Hout;    //Entalpia do fluido na sa�da do elemento de volume [kJ/kg]
        double Houtl;    //Entalpia do fluido na sa�da do elemento de volume levando em conta o incremento de T (dTT) para o c�lculo da derivada [kJ/kg]

        double VazMas;  //vaz�o de fluido pelo elemento [kg/s]
        double VazVol;  //Vaz�o volum�trica [m�/s]
        double vel;     //Velocidade do escoamento como se n�o houvesse o recheio [m/s]
        double Re;      //N�mero de Reinolds para Di�metro da particula e velocidade do vapor no leito sem coque
        double Ra;      //N�mero de Rayleigh
        double Pr;      //N�mero de Prandalt
        double jc;      //Fator de colburn
        double hfrc;    //Coeficiente de troca convec��o em meio poroso [kW/m�.K]
        double hlvr;    //Coeficiente de convec��o livre [kW/m�.K]
        double hli;     //coeficiente de convec��o combindo (livre + for�ada em meio poroso) calculado para (T+dTT) utilizado para derivada
        double hconv;   //Coeficiente de convec��o dentro do c�lculo de ebuli��o [kW/m�.K]
        double hrad;    //Coeficiente de radi��o dentro do c�lculo de ebuli��o [kW/m�.K]
        double hpel;    //Coeficiente de ebuili��o em pel�cula [kW/m�.K]
        double aux;     //Vari�vel auxiliar para o c�lculo do coeficiente de convec��o e ebuli��o em pel�cula
        double n;       //Expoente da regra de mistura do coeficiente de convec��o
        double fTout;   //Fun��o Res�duo da Temperatra de sa�da
        double fToutdT; //Fun��o Res�duo da Temperatra de sa�da levando em conta o incremento de T (dTT) para o c�lculo da derivada
        double flTout;  //Derivada da fun��o Res�duo da Temperatura
        double dh;      //incremento do coeficiente
        double Fh;      //Fun��o Res�duo do coeficiente
        double Fhdh;    //Fun��o Res�duo do coeficiente de sa�da levando em conta o incremento dh
        double Flh;     //Derivada da fun��o Res�duo do coeficiente

        double Qconv;   //Calor calulado pela transferencia de calor por convec��o [kJ/s]
        double QBE;     //Calor calculado pelo balan�o de energia [kJ/s]
        double Qconvl;  //Calor calulado pela transferencia de calor por convec��o levando em conta o incremento de T (dTT) para o c�lculo da derivada [kJ/s]
        double QBEl;    //Calor calculado pelo balan�o de energia levando em conta o incremento de T (dTT) para o c�lculo da derivada [kJ/s]


        //propriedades da �gua no estado f�sico a ser definido
        double roh;
        double cp;
        double K;
        double mu;
        double beta;
        double tens;
        double rohl;    //roh espec�fi�o para l�quidos
        double DHvapl;  //entalpia de vaporiza��o corrigida para corral�o de ebuli��o [kJ/kg]
		
		int aviso;
		int aviso2;
		int aviso3;
		int aviso4;

};

#endif // CONVECCAO_H

//virtual ~CONVECCAO();

    //protected:

    //private:
