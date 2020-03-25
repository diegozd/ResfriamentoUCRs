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
        double h;       //coeficiente de convecção combindo (livre + forçada em meio poroso) [kW/m².K]
        double q;       //Calor de convecção convergido [kJ/s]
		double Qmaximo; //Calor máximo que pode ser trocado pelo fluido.

        string MsgConveccao;

        CONVECCAO();

        double CalcTsai(DADOS&, PROPRIEDADES&, int);
        double calch(double, double, double, double, DADOS&, PROPRIEDADES&, int);
        double calchebu(double, DADOS&, PROPRIEDADES&);
        double calchmax(DADOS&, PROPRIEDADES&);


    private:

        double T0;      //Temperatura de entrada no elemeto finito [K]
        double Tsup;    //Temperatura média da superfício num dado elemeto finito [K]
        double Tout;    //Temperatura de saída do elemnto para convergência [K]
        double TSai;    //Temperatura de saída auxoiliar  do elemento de volume [K]
        double TSail;   //Temepratura de saída + incremento de Temperatura para o cálculo da derivada [K]
        double Tbulkl;  //Temperatura do bulk + incremento de Temperatura para o cálculo da derivada [K]
        double dTT;     //incremento de temperatura [K]
        double DT1;     //variável auxiliar para o cálculo do DTml
        double DT2;     //variável auxiliar para o cálculo do DTml
        double DTml;    //Delta T médio logarítimo [K]
        double DTmll;   //Delta T médio logarítimo levando em conta o incremento de T (dTT) para o cálculo da derivada [K]

        double Hin;     //Entalpia do fluido na estrada do elemento de volume [kJ/kg]
        double Hout;    //Entalpia do fluido na saída do elemento de volume [kJ/kg]
        double Houtl;    //Entalpia do fluido na saída do elemento de volume levando em conta o incremento de T (dTT) para o cálculo da derivada [kJ/kg]

        double VazMas;  //vazão de fluido pelo elemento [kg/s]
        double VazVol;  //Vazão volumétrica [m³/s]
        double vel;     //Velocidade do escoamento como se não houvesse o recheio [m/s]
        double Re;      //Número de Reinolds para Diâmetro da particula e velocidade do vapor no leito sem coque
        double Ra;      //Número de Rayleigh
        double Pr;      //Número de Prandalt
        double jc;      //Fator de colburn
        double hfrc;    //Coeficiente de troca convecção em meio poroso [kW/m².K]
        double hlvr;    //Coeficiente de convecção livre [kW/m².K]
        double hli;     //coeficiente de convecção combindo (livre + forçada em meio poroso) calculado para (T+dTT) utilizado para derivada
        double hconv;   //Coeficiente de convecção dentro do cãlculo de ebulição [kW/m².K]
        double hrad;    //Coeficiente de radição dentro do cãlculo de ebulição [kW/m².K]
        double hpel;    //Coeficiente de ebuilição em película [kW/m².K]
        double aux;     //Variável auxiliar para o cálculo do coeficiente de convecção e ebulição em película
        double n;       //Expoente da regra de mistura do coeficiente de convecção
        double fTout;   //Função Resíduo da Temperatra de saída
        double fToutdT; //Função Resíduo da Temperatra de saída levando em conta o incremento de T (dTT) para o cálculo da derivada
        double flTout;  //Derivada da função Resíduo da Temperatura
        double dh;      //incremento do coeficiente
        double Fh;      //Função Resíduo do coeficiente
        double Fhdh;    //Função Resíduo do coeficiente de saída levando em conta o incremento dh
        double Flh;     //Derivada da função Resíduo do coeficiente

        double Qconv;   //Calor calulado pela transferencia de calor por convecção [kJ/s]
        double QBE;     //Calor calculado pelo balanço de energia [kJ/s]
        double Qconvl;  //Calor calulado pela transferencia de calor por convecção levando em conta o incremento de T (dTT) para o cálculo da derivada [kJ/s]
        double QBEl;    //Calor calculado pelo balanço de energia levando em conta o incremento de T (dTT) para o cálculo da derivada [kJ/s]


        //propriedades da água no estado físico a ser definido
        double roh;
        double cp;
        double K;
        double mu;
        double beta;
        double tens;
        double rohl;    //roh específiço para líquidos
        double DHvapl;  //entalpia de vaporização corrigida para corralão de ebulição [kJ/kg]
		
		int aviso;
		int aviso2;
		int aviso3;
		int aviso4;

};

#endif // CONVECCAO_H

//virtual ~CONVECCAO();

    //protected:

    //private:
