
#include "CONVECCAO.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"
#include "PROPRIEDADES.h"

using namespace std;

CONVECCAO::CONVECCAO()
{
    n = 4.0;      //expoente da regra de mistura do coeficiente de convecção
    dTT = 0.000000001;//dTT = 0.000001;           //incremento de temperatura [K]
    dh = 0.000001;
	aviso = 0;
	aviso2 = 0;
	aviso3 = 0;
	aviso4 = 0;
}

double CONVECCAO::CalcTsai(DADOS &D, PROPRIEDADES &PP, int estado)
{
    T0 = D.Tflu[D.ii][D.jj];
    Tsup = D.Tsup[D.ii][D.jj];
    Hin = PP.Entalpia(T0, D.P, estado);  //'Cálculo da Entalpia de entrada do fluido [kJ/kg]

    if (estado == 1)
    {   //estado vapor
        Tout = D.Tsup[D.ii][D.jj];
        VazMas = D.mvap[D.ii][D.jj];
    }
    else
    {   //estado líquido
        Tout = D.Tsat;
        VazMas = D.mliq[D.ii][D.jj];
    }
	
    h = calch(T0, Tout, Tsup, VazMas, D, PP, estado);
    if (vel > 0)
    {
        TSai = Tsup - exp(-h * D.Apt / (roh * vel * D.Asr * cp)) * (Tsup - T0);  //'estimativa da temperatura na saída do elemento finito segundo 7.93 Incrópera
    }
    else
    {
        TSai = D.Tsat;
    }

    DT2 = Tsup - T0;
    int contloop = 1;       //Variável auxiliar para contagem dos loops
    do
    {
        h = calch(T0, TSai, Tsup, VazMas, D, PP, estado);
        DT1 = Tsup - TSai;
		
        if ((((DT1 > 0) && (DT2 > 0)) || ((DT1 < 0) && (DT2 < 0))) && ((DT1 - DT2) != 0) && (VazMas > 0))
        {
            DTml = (DT1 - DT2) / log(DT1 / DT2);
        }
        else if ((DT1 - DT2) == 0 && (VazMas > 0))
        {
            DTml = Tsup - TSai;
        }
        else
        {
            DTml = 0;
        }

        Qconv = h * D.Apt * DTml;

        if ((TSai > D.Tsat) && (estado == 0))
        {
            Hout = PP.Entalpia(D.Tsat, D.P, estado);
            cp = PP.cpH2O(D.Tsat, estado);
            Hout += VazMas * cp * (TSai - D.Tsat);
        }
        else
        {
            Hout = PP.Entalpia(TSai, D.P, estado);
        }
        QBE = VazMas * (Hout - Hin);
        fTout = Qconv - QBE;
		
        /////////////////////////
        ////    Derivando    ////
        /////////////////////////
        TSail = TSai - dTT;
        Tbulkl = Tbulk - dTT / 2;
        hli = calch(T0, TSail, Tsup, VazMas, D, PP, estado);
        DT1 = Tsup - TSail;
        if ((((DT1 > 0) && (DT2 > 0)) || ((DT1 < 0) && (DT2 < 0))) && ((DT1 - DT2) != 0))
        {
            DTmll = (DT1 - DT2) / log(DT1 / DT2);
        }
        else if ((DT1 - DT2) == 0)
        {
            DTmll = Tsup - TSail;
        }
        else
        {
            DTmll = 0;
        }

        Qconvl = hli * D.Apt * DTmll;
         if ((TSail > D.Tsat) && (estado == 0))
        {
            Houtl = PP.Entalpia(D.Tsat, D.P, estado);
            cp = PP.cpH2O(D.Tsat, estado);
            Houtl += VazMas * cp * (TSail - D.Tsat);
        }
        else
        {
            Houtl = PP.Entalpia(TSail, D.P, estado);
        }
        QBEl = VazMas * (Houtl - Hin);

        fToutdT = Qconvl - QBEl;
        flTout = (fTout - fToutdT) / dTT;

        ///////////////////////////////////////
        ////    Aplicando Newto-Raphson    ////
        ///////////////////////////////////////
        if (flTout != 0)
        {
            Tout = TSai - (fTout / flTout);
        }
        else
        {
            Tout = Tsup;
        }

        ///////////////////////////////////
        ////    Aplicando restições    ////
        ///////////////////////////////////
        if ((Tout - Tsup) > 0)
        {   //Tsaída tem que ser menor ou igual a Tsuperfície
            TSai = Tsup;
            if (((Tout - Tsup) > 0.000001) || (contloop < 3))
            {
                fTout = 1;
            }
        }
        else if((Tout - T0) < 0)
        {   //Tsaída tem que ser maior ou igual a de entrada
            TSai = T0;
            fTout = 1;
        }
        else
        {   //atende as restrições
            TSai = Tout;
        }

        //////////////////////////////////////
        ////    Evitando loop infinito    ////
        //////////////////////////////////////
        if (contloop >= D.nmax)
        {
            MsgConveccao = "Atingiu numero maximo de iteracoes em CalcTsai";
            if ((aviso2 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConveccao << endl;
			aviso2 = 1;
            fTout = 0;
        }
        contloop += 1;
    }while (abs(fTout) > 0.00001);

    q = (Qconv + QBE)/2.0;
	
	Qmaximo = D.Mcoq * D.cpcoq * ((D.Tsup[D.ii][D.jj] + D.Tcoq[D.ii][D.jj])/2.0 - T0) / D.dt;
    
	if (q > Qmaximo) 
	{	
		if ((aviso == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " Calor muito grande retirado do coque." << endl;
		aviso = 1;
		
		q = Qmaximo;
		Hout = q/VazMas + Hin;
		
		fTout = Hout - PP.Entalpia(D.Tsat, D.P, 0);
		fToutdT = Hout - PP.Entalpia(D.Tsat, D.P, 1);
		
		if (((fTout > 0) && (fToutdT < 0)) || ((fTout < 0) && (fToutdT > 0))) 
		{
			TSai = D.Tsat;
		}
		else
		{
			contloop = 1;
			do
			{	
				fTout = Hout - PP.Entalpia(TSai, D.P, estado);
				fToutdT = Hout - PP.Entalpia((TSai + dTT), D.P, estado);
				flTout = (fToutdT-fTout)/dTT;
				TSai = TSai - (fTout / flTout);
				
				if (contloop >= D.nmax)
				{
					MsgConveccao = "Atingiu numero maximo de iteracoes em CalcTsai";
					if ((aviso3 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConveccao << endl;
					aviso3 = 1;
					fTout = 0;
				}
				contloop++;
			}while (abs(fTout) > 0.00001);
		}
	}
	
	Tbulk = (T0 + TSai)/2.0;
    Tfilme = (Tbulk + Tsup)/2.0;

    return TSai;

}

double CONVECCAO::calch(double Tin, double Tout2, double Tsup2, double VazMas2, DADOS &D, PROPRIEDADES &PP, int estado)
{	
    Tbulk = (Tin + Tout2) / 2;
    if ((Tout2 > (4*D.Tsat - 2*Tsup2 - Tin)) && (estado == 0))
    {   //estado == 0 -> liquido
        Tfilme = D.Tsat;
    }
    else
    {
        Tfilme = (Tbulk + Tsup2)/2;
    }
	
    //'Propriedades do vapor
    roh = PP.rohH2O(Tfilme, estado);                                //[kg/m³]
    cp = PP.cpH2O(Tfilme, estado);                                  //[kJ/kg.K]
    K = PP.KH2O(Tfilme, estado);                                    //[kW/m.K]
    mu = PP.muH2O(Tfilme, estado);                                  //[Pa.s]
    beta = PP.betaH2O(Tfilme, estado);                              //[1/K]
	
    VazVol = VazMas2 / roh;                                         //[m³/s]
    vel = VazVol / D.Asr;                                           //[m/s]
    Re = D.Dp * vel * roh / (mu * (1 - D.eps) * D.psi);
    Pr = cp * mu / K;
    if (Re > 0)
    {
        jc = 2.19 * pow(Re,(-2.0/3.0)) + 0.78 * pow(Re,(-0.381));
        hfrc = jc * roh * vel * cp / (pow(Pr,(2.0/3.0)));
        hfrc = hfrc * D.mhc;                                        //[kW/m².K]
    }
    else
    {
        hfrc = 0;                                                   //[kW/m².K]
    }

    if ((Tbulk - Tsup2) > 0)
    {
        hlvr = 0;                                                   //[kW/m².K]
    }
    else
    {
        Ra = (9.81 * beta * (Tsup2 - Tbulk) * pow(D.Dp,3)) / (mu * K / (cp * pow(roh,2)));
        hlvr = (2 + (0.589 * pow(Ra,(1.0/4.0))) / (pow((1 + (0.469 / pow(Pr,(9.0/16.0)))),(4.0/9.0)))) * K / D.Dp;     //[kW/m².K]
        hlvr = hlvr * D.mhl;
    }
	
    aux = pow((pow(hfrc,n) + pow(hlvr,n)),(1.0/n));                 //[kW/m².K]
    return aux;
}

double CONVECCAO::calchebu(double Tsup3, DADOS &D, PROPRIEDADES &PP)
{
    Tfilme = (Tsup3 + D.Tsat)/2;    //[K]
    roh = PP.rohH2O(Tfilme, 1);     //[kg/m³]
    cp = PP.cpH2O(Tfilme, 1);       //[kJ/kg.K]
    K = PP.KH2O(Tfilme, 1);         //[kW/m.K]
    mu = PP.muH2O(Tfilme, 1);       //[Pa.s]
    rohl = PP.rohH2O(D.Tsat, 0);    //[kg/m³]

    DHvapl = PP.DHVapTsat + 0.8 * cp * (Tsup3 - D.Tsat);        //Correção da entalpia de vaporização [kJ/kg]
    hconv = (0.67 * K / D.Dp) * pow(((roh * 9.81 * (rohl - roh) * DHvapl * pow(D.Dp,3)) / (mu * K * (Tsup3 + D.Tsat))),(1.0/4.0));    //'h ebulição em película Eq. 10.9 Incropera [kW/m²K]

    hrad = D.emiss * 5.67 * pow(10 ,-8) * (pow(Tsup3,4) - pow(D.Tsat,4)) / (Tsup3 - D.Tsat) / 1000;      //h de radiação [kW/m²K]

    hpel = hconv + (3.0/4.0) * hrad;
    int contloop = 1;       //Variável auxiliar para contagem dos loops
    do
    {
        Fh = (pow(hconv,(4.0/3.0)) + hrad * pow(hpel,(4.0/3.0))) - pow(hpel,(4.0/3.0));
        Fhdh = (pow(hconv,(4.0/3.0)) + hrad * pow((hpel + dh),(4.0/3.0))) - pow((hpel + dh),(4.0/3.0));
        Flh = (Fhdh - Fh) / dh;
        hpel = hpel - Fh / Flh;

        //'evitando loop infinito
        if (contloop >= D.nmax)
        {
            MsgConveccao = "Atingiu numero maximo de iteracoes em calchebu";
            if ((aviso4 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConveccao << endl;
			aviso4 = 1;
            Fh = 0;
        }
        contloop += 1;
    }while (abs(Fh) > 0.00000001);

    aux = hpel * D.mhe;
    return aux;
}

double CONVECCAO::calchmax(DADOS &D, PROPRIEDADES &PP)
{
    roh = PP.rohH2O(D.Tsat, 1);     //[kg/m³]
    rohl = PP.rohH2O(D.Tsat, 0);    //[kg/m³]
    tens = PP.sigH2O(D.Tsat);

    aux = (0.149 * PP.DHVapTsat * roh * pow((tens * 9.81 * (rohl - roh) / pow(roh,2)),(1.0/4.0))) / 30;
    return aux;
}
