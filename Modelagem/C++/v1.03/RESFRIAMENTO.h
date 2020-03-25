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

        double TotalAgua;                   //totalizador da água que foi injetada no leito [kg]
        double VaporGeradoTotal;            //totalizador de todo vapor gerado no leito [kg]
        double AgAcumuAcimaLeito;           //totalizador da água acumulada acima do leito [kg]
        double AgAcumu;                     //totalizador da água acumulada no leito [kg]
        double PicoVazao;                   //Pico de vazão de vapor [kg/s]

        vector<vector<double> > GeraVap;     //Geração de vapor no elemento de volume [kg/s]
        vector<double> NIVEL;               //Nível de líquido no leito ao longo do tempo [m]
        vector<double> DistEspuma;          // vetor de distribuição da espuma

        string MsgResfriamento;

        RESFRIAMENTO(DADOS&, PROPRIEDADES&);

        void espuma(DADOS&);
        double mistura_correntes(double, double, double, double, DADOS&, PROPRIEDADES&, int);

    private:

        double Tsai;        //Temperatura de saída auxiliar do elemento de volume [K]
        double DT1;         //variável auxiliar para o cálculo do DTml [K]
        double DT2;         //variável auxiliar para o cálculo do DTml [K]
        double DTml;        //Delta T médio logarítimo [K]
        double DTe;         //Delta Te para ebulição[K]
        double Hsat;        //Entalpia do fluido líquido na temperatura de saturação [kJ/kg]
        double Qsensivel;   //Calor sensível [kJ/s]
        double Qebu;        //Calor de ebulição [kJ/s]
        double Qconv;       //Calor de convgecção mista [kJ/s]
        double Qcoq;        //Calor de ebulição retirado do coque [kJ/s]
        double hebu;        //Coeficiente de ebuilição em película [kW/m².K]
        double hconv;       //Coeficiente de convecção [kW/m².K]
        double hmin;        //Coeficiente mínimo de ebuilição em película [kW/m².K]
        double a1;          //coeficiente angular da reta de hebu para ebulição nucleada e transição
        double b1;          //coeficiente linear da reta de hebu para ebulição nucleada e transição
        double rohl;        //roh específiço para líquidos [kg/m³]
        double rohv;        //roh específiço para vapor [kg/m³]
        double Tvap;        //Temperatura auxiliar do vapor [K]
        double DVzOld;      //VAriável auxiliar
        double Dmed;        //Diâmetro médio [m]
        double VolCoq;      //'Volume do coque no elemento [m³]
        double L;           //variável de acrescimo do nível [m]
        double L2;          //variável auxioliar de acrescimo do nível [m]
        double Lvelho;      //variável auxiliar para o cálculo do acrescimo do nível [m]
        double DIFF2;       //variável de convergência do nivel
        double ContaCurva;  //contador da curva de resfriamento
        double HliqIn;      //entalpia do líquido na entrada do elemento
        double DeslocVol;   //vazão deslocada pelo líquido não vaporizado [kg/s]
        double Dd;          //Delta de diâmetro no cone [m]
        double hmax;        //Coeficiente máximo de ebuilição em película [kW/m².K]
        double VazLiqEl;    //auxiliar para vazão de líquido no elemento [kg/s]

        int n;              //variável auxiliar para espuma
        int aux;            //variável auxliar para espuma 2
		
		int aviso;
		int aviso2;


};

#endif // RESFRIAMENTO_H

//virtual ~RESFRIAMENTO();

    //protected:

    //private:
