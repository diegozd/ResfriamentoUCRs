
#include "PROPRIEDADES.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"

using namespace std;
PROPRIEDADES::PROPRIEDADES(DADOS &D)
{
    //Calcula a Tempertuara de satura��o da s�gua em K dado P em bar abs
    Tsat = 1 / (0.00895681158996 + -0.00627384806376 *  pow(D.P,0.0318661188564));
    DHVapTsat = DHvapH2O(Tsat);
    D.Tsat = Tsat;
	
	aviso_PSO = D.PSO;
	aviso_roh = 0;
	aviso_cp = 0;
	aviso_K = 0;
	aviso_mu = 0;
	aviso_beta = 0;
	aviso_Entalp = 0;
};

//Retorna massa espec�fica da �gua em kg/m� nas fases l�quida ou vapor
double PROPRIEDADES::rohH2O(double T, int estado)
{
    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }


    if (vai == 0)       //liquido
    {
        //massa espec�fica da �gua saturada @ T  [kg/m�]
        aux = -0.002452231 * pow(T,2) + 1.124941 * T + 879.9971;
		if(T < 273.15) aux = -0.002452231 * pow(273.15,2) + 1.124941 * 273.15 + 879.9971;
    }
    else if (vai == 1)  //vapor
    {
        //calcula massa espec�fica do vapor a 2barg em kg/m�
        aux = 0.000004166821 * pow(T,2) - 0.007054861 * T + 3.821986;
    }
    else                //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em rohH2O";
        if ((aviso_roh == 0) && (aviso_PSO == 0)) cout << MsgPropriedades << endl;
		aviso_roh = 1;
    }

    return aux;

};


//Retorna capacidade calor�fica da �gua em kJ/kg.K nas fases l�quida ou vapor
double PROPRIEDADES::cpH2O(double T, int estado)
{
    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }


    if (vai == 0)           //liquido
    {
        //capacidade calor�fica da �gua saturada @ T [kJ/kg.K]
        aux = 0.00001178632 * pow(T,2) - 0.007445843 * T + 5.354467;
		if (T < 273.15) aux = 0.00001178632 * pow(273.15,2) - 0.007445843 * 273.15 + 5.354467;
    }
    else if (vai == 1)      //vapor
    {
        //calcula capacidade calor�fica do vapor a 1 atm [kJ/kg.K]
        aux = 6.62561202404 -0.0290862268271 * T + 6.61776463052E-05 * pow(T,2) + -6.51294807743E-08 * pow(T,3) + 2.38881298235E-11 * pow(T,4);
    }
    else                    //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em cpH2O";
        if ((aviso_cp == 0) && (aviso_PSO == 0)) cout << MsgPropriedades << endl;
		aviso_cp = 1;
    }

    return aux;

};


//Retorna condutividade t�rmicado da �gua em kW/m.K nas fases l�quida ou vapor
double PROPRIEDADES::KH2O(double T, int estado)
{
    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }


    if (vai == 0)           //liquido
    {
        //condutividade t�rmicado da �gua saturada @ T [kW/m.K]
        aux = -0.000000006532962 * pow(T,2) + 0.000005318497 * T - 0.0003936951;
		if (T < 273.15) aux = -0.000000006532962 * pow(273.15,2) + 0.000005318497 * 273.15 - 0.0003936951;
    }
    else if (vai == 1)      //vapor
    {
        //calcula capacidade condutividade t�rmicado vapor a 2barg em kW/m.K
        aux = 0.0000001088875 * T - 0.00001781313;
    }
    else                    //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em KH2O";
        if ((aviso_K == 0) && (aviso_PSO == 0)) cout << MsgPropriedades << endl;
		aviso_K = 1;
    }

    return aux;

};


//Retorna viscosidade da �gua em Pa.s nas fases l�quida ou vapor
double PROPRIEDADES::muH2O(double T, int estado)
{
    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }


    if (vai == 0)           //liquido
    {
        //viscosidade da �gua saturada @ T [Pa.s]
        aux = 1 / (-8883.00826068 + 45.9874310456 * pow(T,0.945111886617));
		if (T < 273.15) aux = 1 / (-8883.00826068 + 45.9874310456 * pow(273.15,0.945111886617)); // restri��o da viscosidade
    }
    else if (vai == 1)      //vapor
    {
        //calcula capacidade viscosidade vapor a 2barg em Pa.s == kg/ms
        aux = 0.00000004140873 * T - 0.000003488596;
    }
    else                    //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em muH2O";
        if ((aviso_mu == 0) && (aviso_PSO == 0)) cout << MsgPropriedades << endl;
		aviso_mu = 1;
    }

    return aux;

};


//Retorna tens�o superfuicial da �gua l�quida em N/m
double PROPRIEDADES::sigH2O(double T)
{
    // tens�o superfuicial da �gua saturada @ T [N/m]
    aux = -0.0000002057757 * pow(T,2) - 0.0000380553 * T + 0.1017064;

    return aux;

};


//Retorna coeficiente de expans�o da �gua em 1/K nas fases l�quida ou vapor
double PROPRIEDADES::betaH2O(double T, int estado)
{
    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }


    if (vai == 0)           //liquido
    {
        //coeficiente de expans�o da �gua (beta) saturada @ T [1/K]
        aux = -0.00000001442629 * pow(T,2) + 0.00001631945 * T - 0.003320243;
		if (T < 273.15) aux = -0.00000001442629 * pow(273.15,2) + 0.00001631945 * 273.15 - 0.003320243; // restri��o de temperatura
    }
    else if (vai == 1)      //vapor
    {
        //coeficiente de expans�o da vapor como se fosse gas idea @ T [1/K]
        aux = 1/T;
    }
    else                    //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em betaH2O";
        if ((aviso_beta == 0) && (aviso_PSO == 0)) cout << MsgPropriedades << endl;
		aviso_beta = 1;
    }

    return aux;

};


//Retorna entalpia de vaporiza��o da �gua em kJ/kg
double PROPRIEDADES::DHvapH2O(double T)
{
    //Calcula a entalpia de vaporiza��o da �gua em kJ/kg dado T em K
    aux = 2787.26101875 -0.126009015962 * T + -0.0034734972292 * pow(T,2);
	if (T < 273.15) aux = 2787.26101875 -0.126009015962 * 273.15 + -0.0034734972292 * pow(273.15,2);

    return aux;

};


//Retorna entalpia da �gua em kJ/kg
double PROPRIEDADES::Entalpia(double T, double P, int estado)
{
    //limites de T
    double Tref = 873.15;       //Temperatura de refer�ncia que � o limite superior da correla��o do cp [K]
    double Tlim = 380;          //limite inferior da correla��o do cp [K]

    //contantes do Cp vapor
    long double A = 6.62561202404;
    long double B = -0.0290862268271;
    long double C = 6.61776463052E-05;
    long double D = -6.51294807743E-08;
    long double E = 2.38881298235E-11;

    //constantes Cp liquido
    long double Al = 0.00001178632;
    long double Bl = -0.007445843;
    long double Cl = 5.354467;

    //Constantes para estimas a varia��o da entalpia do g�s com a press�o
    long double a1 = -6.10703633699;
    long double b1 = 2174.88990338;
    long double c1 = 0.506192451569;
    long double a2 = -5.97456441139;
    long double b2 = 2178.87567477;
    long double c2 = 0.488011821528;

    //auxiliares
    double Hatm;
    double Ap;
    double Bp;
    double DHp;
    double DHvap;
    double Hliq;

    if (T < Tsat)
    {
        vai = 0; //l�quido
    }
    else if (T > Tsat)
    {
        vai = 1; //vapor
    }
    else
    {
        vai = estado; //bif�sico
    }

    if (vai == 0)           //liquido
    {
        if (T < Tlim) //'limite inferior da correla��o de cp do vapor
        {
            Hatm = A * (Tlim - Tref) + B / 2 * (pow(Tlim,2) - pow(Tref,2)) + C / 3 * (pow(Tlim,3) - pow(Tref,3)) + D / 4 * (pow(Tlim,4) - pow(Tref,4)) + E / 5 * (pow(Tlim,5) - pow(Tref,5)); //Entapia a 1amt [kJ/kg]

            Ap = exp(a1 + b1 / Tlim + c1 * log(Tlim));  // 'coeficiente angular da reta deltaH vs P
            Bp = exp(a2 + b2 / Tlim + c2 * log(Tlim));  // 'coeficiente linear da reta deltaH vs P
            DHp = Ap * P - Bp;                           //'corre��od a entalpia com a press�o

            DHvap = DHvapH2O(Tlim);

            Hliq = Al / 3 * (pow(Tlim,3) - pow(T,3)) + Bl / 2 * (pow(Tlim,2) - pow(T,2)) + Cl * (Tlim - T);

            aux = Hatm - DHp - DHvap - Hliq;        //Entapila [kJ/kg]
        }
        else
        {
            Hatm = A * (T - Tref) + B / 2 * (pow(T,2) - pow(Tref,2)) + C / 3 * (pow(T,3) - pow(Tref,3)) + D / 4 * (pow(T,4) - pow(Tref,4)) + E / 5 * (pow(T,5) - pow(Tref,5)); //'Entapia a 1amt [kJ/kg]

            Ap = exp(a1 + b1 / T + c1 * log(T)); //'coeficiente angular da reta deltaH vs P
            Bp = exp(a2 + b2 / T + c2 * log(T)); //'coeficiente linear da reta deltaH vs P
            DHp = Ap * P - Bp;                   //'corre��od a entalpia com a press�o [kJ/kg]

            DHvap = DHvapH2O(T);                 //'Entalpia de vaporiza��o a T [kJ/kg]

            aux = Hatm - DHp - DHvap;            //'Entapila [kJ/kg]
        }
    }
    else if (vai == 1)      //vapor
    {
        Hatm = A * (T - Tref) + B / 2 * (pow(T,2) - pow(Tref,2)) + C / 3 * (pow(T,3) - pow(Tref,3)) + D / 4 * (pow(T,4) - pow(Tref,4)) + E / 5 * (pow(T,5) - pow(Tref,5)); //'Entapia a 1amt [kJ/kg]

        Ap = exp(a1 + b1 / T + c1 * log(T)); //'coeficiente angular da reta deltaH vs P
        Bp = exp(a2 + b2 / T + c2 * log(T)); //'coeficiente linear da reta deltaH vs P
        DHp = Ap * P - Bp;                   //'corre��od a entalpia com a press�o[kJ/kg]

        aux = Hatm - DHp;                    //'Entapila [kJ/kg]
    }
    else                    //erro
    {
        aux = 0;
        MsgPropriedades = "Erro em Entalpia";
        if ((aviso_Entalp == 0) && (aviso_PSO == 0)) cout << " " << MsgPropriedades << endl;
		aviso_Entalp = 1;
    }

    return aux;

}
