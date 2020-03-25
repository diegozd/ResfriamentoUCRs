
#include "CONDUCAO.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "CONVECCAO.h"

using namespace std;

CONDUCAO::CONDUCAO()
{
    dksi = 0.000000000001;
    dttt = 0.000000000001;
	//aviso = 0;
	aviso2 = 0;
	aviso3 = 0;
	aviso4 = 0;
    //ctor
}

double CONDUCAO::CalcTsup(DADOS &D, PROPRIEDADES &PP, CONVECCAO &HLV)
{
    Tcoq0 = D.Tcoq[D.ii][D.jj];
    Bi = HLV.h * (D.Dp/2) / D.Kcoq;
    alphaCoq = D.Kcoq / (D.rohcoq * D.cpcoq);
    Fo = alphaCoq * D.dt / (pow((D.Dp/2),2));

    n = 1;
    tetaold = 0;
    int contloop = 1;       //Variável auxiliar para contagem dos loops
    do
    {
        ksi = raizksi(Bi, n, D);
        Cn = 4 * (sin(ksi) - ksi * cos(ksi)) / (2 * ksi - sin(2 * ksi));
        teta = tetaold + (Cn * exp(-Fo * pow(ksi,2)) * (sin(ksi) / ksi));
		
        DIFF1 = abs((teta - tetaold) / teta);
        tetaold = teta;
        n += 1;

        //////////////////////////////////////
        ////    Evitando loop infinito    ////
        //////////////////////////////////////
        if (contloop == D.nmax * 10)
        {
            MsgConducao = "Atingiu numero maximo de iteracoes em CalcTsup";
            if ((aviso2 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConducao << endl;
			aviso2 = 1;
            DIFF1 = 0;
        }
        contloop += 1;

    }while (DIFF1 > 0.00001);

    Tsuper = teta * (Tcoq0 - HLV.Tbulk) + HLV.Tbulk;

    return Tsuper;

}

double CONDUCAO::CalcTcoq(DADOS &D, PROPRIEDADES &PP, CONVECCAO &HLV)
{
    Tsup0 = D.Tsup[D.ii][D.jj];
    QQ = HLV.q * D.dt;
    DTcoq = QQ / (D.Mcoq * D.cpcoq);

    Tcok = Tcoq0 - dttt;
    int contloop = 1;       //Variável auxiliar para contagem dos loops
    do
    {
        DT1 = Tsup0 - Tsuper;
        DT2 = Tcoq0 - Tcok;

		if ((((DT1 > 0) && (DT2 > 0)) || ((DT1 < 0) && (DT2 < 0))) && (abs(DT1 - DT2) >= 0.0000000001))
        {
            fT = DTcoq - (DT1 - DT2) / log(DT1 / DT2);
            fTdT = DTcoq - (DT1 - (Tcoq0 - (Tcok - dttt))) / log(DT1 / (Tcoq0 - (Tcok - dttt)));
            flT = (fT - fTdT) / dttt;
            Tcok = Tcok - fT / flT;
        }
        else if (abs(DT1 - DT2) < 0.0000000001)
        {
            fT = DTcoq - DT2;
            fTdT = DTcoq - DT2 - dttt;
            flT = (fT - fTdT) / dttt;
            Tcok = Tcok - fT / flT;
        }
        else
        {
            Tcok = Tcoq0;
            fT = 0;
        }
		
		/*if ((Tcok < Tsuper) && (contloop > 2) && (D.ii != 0))
		{
			if ((aviso == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << "Tcoq menor que Tsuperficie." << endl;
			aviso = 1;
			
			Tcok = Tsuper;
			
			if (contloop2 > 3) fT = 0;
			
			contloop2 ++;
		}*/

        //////////////////////////////////////
        ////    Evitando loop infinito    ////
        //////////////////////////////////////
        if (contloop == D.nmax)
        {
            MsgConducao = "Atingiu numero maximo de iteracoes em CalcTcoq";
            if ((aviso3 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConducao << endl;
			aviso3 = 1;
            fT = 0;
        }
        contloop += 1;
    }while (abs(fT) > 0.0001); //0.001 //0.000001
	
	if (((Tcok - Tsup0) < 0 ) && D.ii != 0 )
	{
		Tcok = Tsup0;
	}

    return Tcok;

}



double CONDUCAO::raizksi(double Bi,int n, DADOS &D)
{
    ksi = (n - 1)*D.pi + D.pi/2;
    int contloop = 1;       //Variável auxiliar para contagem dos loops
    do
    {
        fksi = 1 - ksi * cos(ksi) / sin(ksi) - Bi;

        /////////////////////////
        ////    Derivando    ////
        /////////////////////////
        ksidksi = ksi - dksi;
        fksidksi = 1 - ksidksi * cos(ksidksi) / sin(ksidksi) - Bi;
        flksi = (fksi - fksidksi) / dksi;

        if (flksi == 0)
        {
            if (contloop == 1)
            {
                ksi = (n - 1) * D.pi + D.pi;
                fksi = 1 - ksi * cos(ksi) / sin(ksi) - Bi;
                ksidksi = ksi - dksi;
                fksidksi = 1 - ksidksi * cos(ksidksi) / sin(ksidksi) - Bi;
                flksi = (fksi - fksidksi) / dksi;
            }
            else
            {
                MsgConducao = "Erro em raizksi";
                cout << D.ii << " " << D.jj << " " << MsgConducao << endl;
            }
        }

        ///////////////////////////////////////
        ////    Aplicando Newto-Raphson    ////
        ///////////////////////////////////////
        aux = ksi - fksi / flksi;

        ///////////////////////////////////
        ////    Aplicando restições    ////
        ///////////////////////////////////
        if (aux > (n * D.pi))
        {
            ksi = n * D.pi;
        }
        else
        {
            ksi = aux;
        }

        //////////////////////////////////////
        ////    Evitando loop infinito    ////
        //////////////////////////////////////
        if (contloop == D.nmax)
        {
            MsgConducao = "Atingiu numero maximo de iteracoes em raizksi";
            if ((aviso4 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgConducao << endl;
			aviso4 = 1;
            fksi = 0;
        }
        contloop += 1;

    }while (abs(fksi) > 0.0001); //0.0001 //0.01

    return aux;

}
