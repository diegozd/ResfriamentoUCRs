#ifndef RESFRIAMENTO_H
#define RESFRIAMENTO_H

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "CONVECCAO.h"
#include "CONDUCAO.h"

using namespace std;

class RESFRIAMENTO
{
    public:

        int nEsp;                           //numero de elementos de volume qua a espuma abrange

        double TotalAgua;                   //totalizador da �gua que foi injetada no leito [kg]
        double VaporGeradoTotal;            //totalizador de todo vapor gerado no leito [kg]
        double AgAcumuAcimaLeito;           //totalizador da �gua acumulada acima do leito [kg]
        double AgAcumu;                     //totalizador da �gua acumulada no leito [kg]
        double PicoVazao;                   //Pico de vaz�o de vapor [kg/s]

        vector<vector<double> > GeraVap;     //Gera��o de vapor no elemento de volume [kg/s]
        vector<double> NIVEL;               //N�vel de l�quido no leito ao longo do tempo [m]
        vector<double> DistEspuma;          // vetor de distribui��o da espuma

        string MsgResfriamento;

        RESFRIAMENTO(DADOS&, PROPRIEDADES&);

        void espuma(DADOS&);
        double mistura_correntes(double, double, double, double, DADOS&, PROPRIEDADES&, int);

    private:

        double Tsai;        //Temperatura de sa�da auxiliar do elemento de volume [K]
        double DT1;         //vari�vel auxiliar para o c�lculo do DTml [K]
        double DT2;         //vari�vel auxiliar para o c�lculo do DTml [K]
        double DTml;        //Delta T m�dio logar�timo [K]
        double DTe;         //Delta Te para ebuli��o[K]
        double Hsat;        //Entalpia do fluido l�quido na temperatura de satura��o [kJ/kg]
        double Qsensivel;   //Calor sens�vel [kJ/s]
        double Qebu;        //Calor de ebuli��o [kJ/s]
        double Qconv;       //Calor de convgec��o mista [kJ/s]
        double Qcoq;        //Calor de ebuli��o retirado do coque [kJ/s]
        double hebu;        //Coeficiente de ebuili��o em pel�cula [kW/m�.K]
        double hconv;       //Coeficiente de convec��o [kW/m�.K]
        double hmin;        //Coeficiente m�nimo de ebuili��o em pel�cula [kW/m�.K]
        double a1;          //coeficiente angular da reta de hebu para ebuli��o nucleada e transi��o
        double b1;          //coeficiente linear da reta de hebu para ebuli��o nucleada e transi��o
        double rohl;        //roh espec�fi�o para l�quidos [kg/m�]
        double rohv;        //roh espec�fi�o para vapor [kg/m�]
        double Tvap;        //Temperatura auxiliar do vapor [K]
        double DVzOld;      //VAri�vel auxiliar
        double Dmed;        //Di�metro m�dio [m]
        double VolCoq;      //'Volume do coque no elemento [m�]
        double L;           //vari�vel de acrescimo do n�vel [m]
        double L2;          //vari�vel auxioliar de acrescimo do n�vel [m]
        double Lvelho;      //vari�vel auxiliar para o c�lculo do acrescimo do n�vel [m]
        double DIFF2;       //vari�vel de converg�ncia do nivel
        double ContaCurva;  //contador da curva de resfriamento
        double HliqIn;      //entalpia do l�quido na entrada do elemento
        double DeslocVol;   //vaz�o deslocada pelo l�quido n�o vaporizado [kg/s]
        double Dd;          //Delta de di�metro no cone [m]
        double hmax;        //Coeficiente m�ximo de ebuili��o em pel�cula [kW/m�.K]
        double VazLiqEl;    //auxiliar para vaz�o de l�quido no elemento [kg/s]

        int n;              //vari�vel auxiliar para espuma
        int aux;            //vari�vel auxliar para espuma 2
		
		int aviso;
		int aviso2;


};

#endif // RESFRIAMENTO_H

//virtual ~RESFRIAMENTO();

    //protected:

    //private:
