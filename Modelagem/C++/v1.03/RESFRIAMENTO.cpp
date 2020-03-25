//Para usar o comando | system("PAUSE"); | é necessário inluir os seguintes comandos
#include <stdlib.h>

#include "RESFRIAMENTO.h"

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

RESFRIAMENTO::RESFRIAMENTO(DADOS &D, PROPRIEDADES &PP)
{
	aviso = 0;
	aviso2 = 0;
	double auxNIVEL = 0;
	vector <int> CalculeiInterface;
	
    //criando as classes de convecção e condução para líquido e vapor
    CONVECCAO HLV;
    CONDUCAO KLV;

    //Calculando constantes
    D.NumL = D.AltLeito / D.dL;                                  //numero de elementos de espaço
    D.NumT = D.CurvResf[(D.CurvResf.size()-1)][1] * 60 / D.dt;   //numero de elementos no tempo
    Dd = (D.Dc - D.Dcone) / (D.hcone / D.dL);                  //Delta de diâmetro no cone [m]
    hmax = HLV.calchmax(D, PP);                                //Coeficiente máximo de ebuilição em película [kW/m².K]
    
	//ajustando o numero de intervalos de espaço 
	if ((D.AltLeito / D.dL) > D.NumL ) { D.NumL = D.NumL + 1; }
	
    //atribuindo valores iniciais a algumas variáveis
    Lvelho = 0;  //variável auxiliar para o cálculo do acrescimo do nível

    //dimensionando matrizes
    D.Tflu.resize(D.NumT);
    D.Tcoq.resize(D.NumT+1);
    D.Tsup.resize(D.NumT+1);
    D.mliq.resize(D.NumT);
    D.mvap.resize(D.NumT);
    GeraVap.resize(D.NumT);
    NIVEL.resize(D.NumT+1);
	CalculeiInterface.resize(D.NumT);

    D.Tcoq[0].resize(D.NumL);
    D.Tsup[0].resize(D.NumL);
    for (int i = 0 ; i < D.NumT ; i++)
    {
        D.Tflu[i].resize(D.NumL+1);
        D.Tcoq[i+1].resize(D.NumL);
        D.Tsup[i+1].resize(D.NumL);
        D.mliq[i].resize(D.NumL+1);
        D.mvap[i].resize(D.NumL+1);
        GeraVap[i].resize(D.NumL);
    }

	//										//
	///									   ///
    ////                                  ////
    /////                                /////
    //////  INÍCIO LÓGICO DO PROGRAMA   //////
    /////                                /////
    ////                                  ////
	///									   ///
	//										//
	
    NIVEL[0] = 0;
    TotalAgua = 0;
    VaporGeradoTotal = 0;
    AgAcumuAcimaLeito = 0;
    AgAcumu = 0;
    PicoVazao = 0;

    espuma(D);

    for (D.ii = 0 ; D.ii < D.NumT ; D.ii++)   //NumT
    {  //início do loop de tempo
		
		CalculeiInterface[D.ii] = 0;
		
        ContaCurva = 0;
        D.Tflu[D.ii][0] = 0.0;
        while (D.Tflu[D.ii][0] == 0.0)
        {   //Lendo vazão e temperatura de entrada no leito - roda até conseguir ler a primeira temperatura
            if (D.ii*D.dt < D.CurvResf[ContaCurva][1]*60)
            {   //ler valores se o Tempo for adequado
                D.Tflu[D.ii][0] = D.CurvResf[ContaCurva][2] + 273.15;           //temperatura do fluido na entrada do elemento de volume [K]
                if (D.Tflu[D.ii][0] <= PP.Tsat)
                {   //entrada líquida
                    D.mliq[D.ii][0] = D.CurvResf[ContaCurva][0] / 3.6;  //vazão mássica de líquido na entrada do elemento de volume [kg/s]
                    D.mvap[D.ii][0] = 0.0;                                      //vazão mássica de gás na entrada do elemento de volume [kg/s]
                }
                else
                {   //entrada vapor
                    D.mliq[D.ii][0] = 0.0;                                      //vazão mássica de líquido na entrada do elemento de volume [kg/s]
                    D.mvap[D.ii][0] = D.CurvResf[ContaCurva][0]  / 3.6;  //vazão mássica de gás na entrada do elemento de volume [kg/s]
                }
            }

            ContaCurva += 1;
        }
		
        n = 0;
        TotalAgua +=  D.mliq[D.ii][0] * D.dt;       //totalizador de água injetada no leito [kg]

        for(D.jj = 0 ; D.jj < D.NumL ; D.jj++)
        {   //início do loop de espaço

            //inicializando temperaturas no interior e na superfície do coque no tempo 0
            if (D.ii == 0)
            {
                D.Tcoq[0][D.jj] = D.ToCoq;
                D.Tsup[0][D.jj] = D.ToCoq;
            }

            //cálculo do diâmetro médio no elemento de volume
            if(D.jj*D.dL < D.hcone)
            {
                Dmed = D.Dcone + D.jj*Dd + Dd/2;      //Diâmetro médio do cone [m]
            }
            else
            {
                Dmed = D.Dc;
            }

            D.Asr = D.pi * pow(Dmed,2) / 4;
            D.Apt = 3 * (1 - D.eps) * D.pi * D.dL * pow((Dmed/2),2) / (D.Dp/2);
            VolCoq = D.dL * D.pi * pow(Dmed,2) / 4 * (1 - D.eps);
            D.Mcoq = VolCoq * D.rohcoq;

            if(n == 0)
            {
                aux = D.jj;
            }

            ///////////////////////////////////////////////
            ////    ELEMENTO COMPLETAMENTE SUBMERSO    ////
            ///////////////////////////////////////////////
            if(NIVEL[D.ii] >= (D.jj*D.dL + D.dL))
            {	
                HliqIn = PP.Entalpia(D.Tflu[D.ii][D.jj],D.P,0);       //'Entalpia do líquido na entrada do elemento de volume [kJ/kg]
                Tsai = HLV.CalcTsai(D, PP, 0); 
				
                if (Tsai >= D.Tsat)
                {   //houve ebulição
                    D.Tflu[D.ii][D.jj+1] = D.Tsat;
					
                    DT1 = D.Tsup[D.ii][D.jj] - D.Tflu[D.ii][D.jj];
                    DT2 = D.Tsup[D.ii][D.jj] - D.Tflu[D.ii][D.jj+1];
                    if ((((DT1 > 0) && (DT2 > 0)) || ((DT1 < 0) && (DT2 < 0))) && ((DT1 - DT2) != 0))
                    {
                        DTml = (DT1 - DT2) / log(DT1 / DT2);
                    }
                    else if ((DT1 - DT2) == 0)
                    {
                        DTml = D.Tsup[D.ii][D.jj] - D.Tsat;
                    }
                    else
                    {
                        DTml = 0;
                    }

                    Hsat = PP.Entalpia(D.Tsat, D.P, 0);
                    Qsensivel = D.mliq[D.ii][D.jj] * (Hsat - HliqIn);
                    DTe = D.Tsup[D.ii][D.jj] - D.Tsat;

                    if (DTe >= 120)
                    {   //ebulição em película
                        hebu = HLV.calchebu(D.Tsup[D.ii][D.jj], D, PP);
                        Qebu = hebu * D.Apt * DTe;
                        Qconv = Qsensivel;
                    }
                    else if (DTe >= 30)
                    {   //ebulição em regime de transição
                        hmin = HLV.calchebu((D.Tsat+120), D, PP);
                        a1 = (hmax - hmin) / (30 - 120);
                        b1 = hmax - a1 * 30;

                        hebu = (a1 * DTe + b1) * D.mhe;
                        Qebu = hebu * D.Apt * DTe;
                        Qconv = Qsensivel;
                    }
                    else if (DTe >= 5)
                    {   //ebulição nucleada
                        hmin = HLV.calch(D.Tflu[D.ii][D.jj], D.Tsat, (D.Tsat + 5), D.mliq[D.ii][D.jj], D, PP, 0);
                        a1 = (hmax - hmin) / (30 - 5);
                        b1 = hmax - a1 * 30;

                        hebu = (a1 * DTe + b1) * D.mhe;
                        Qebu = hebu * D.Apt * DTe;
                        Qconv = Qsensivel;
                    }
                    else
                    {   //ebulição livre
                        Qebu = 0;
                        Qconv = (HLV.h * D.Apt * DTml);
                    }

                    Qcoq = Qconv + Qebu;
                    hconv = Qcoq / (D.Apt * DTml);
                    GeraVap[D.ii][D.jj] = (Qcoq - Qsensivel) / PP.DHVapTsat;
                    if (GeraVap[D.ii][D.jj] > D.mliq[D.ii][D.jj])
                    {
                        GeraVap[D.ii][D.jj] = D.mliq[D.ii][D.jj];
                        Qebu = D.mliq[D.ii][D.jj] * PP.DHVapTsat;
                        Qcoq = Qconv + Qebu;
                        hconv = Qcoq / (D.Apt * DTml);
                    }
                }
                else
                {   //não houve ebulição
                    D.Tflu[D.ii][D.jj+1] = Tsai;
                    hconv = HLV.h;
                    Qcoq = HLV.q;
                    GeraVap[D.ii][D.jj] = 0;
                }

                D.mliq[D.ii][D.jj+1] = D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj];
                D.mvap[D.ii][D.jj+1] = D.mvap[D.ii][D.jj] + GeraVap[D.ii][D.jj];

                HLV.h = hconv;
                HLV.Tbulk = (D.Tflu[D.ii][D.jj] + D.Tflu[D.ii][D.jj+1])/2;
                HLV.q = Qcoq;

                D.Tsup[D.ii+1][D.jj] = KLV.CalcTsup(D, PP, HLV);
                D.Tcoq[D.ii+1][D.jj] = KLV.CalcTcoq(D, PP, HLV);
				
				/*if (D.ii == 4201 && D.jj == 59)
				{
					cout << "deu ruim - 1" << endl;
					cout << "D.ii = " << D.ii << " D.jj = " << D.jj << endl;
					cout << "GeraVap[D.ii][D.jj] = " << GeraVap[D.ii][D.jj] << endl;
					cout << "D.mliq[D.ii][D.jj+1] = " << D.mliq[D.ii][D.jj+1] << endl;
					cout << "D.mvap[D.ii][D.jj+1] = " << D.mvap[D.ii][D.jj+1] << endl;
					cout << "D.Tflu[D.ii][D.jj+1] = " << D.Tflu[D.ii][D.jj+1] << endl;
					cout << "D.Tsup[D.ii+1][D.jj] = " << D.Tsup[D.ii+1][D.jj] << endl;
					cout << "D.Tcoq[D.ii+1][D.jj] = " << D.Tcoq[D.ii+1][D.jj] << endl;
					
					system("PAUSE");
				}*/

                if (GeraVap[D.ii][D.jj] > 0)
                {
					auxNIVEL = NIVEL[D.ii];
					
                    rohl = PP.rohH2O(D.Tflu[D.ii][D.jj], 0);
                    //AgAcumu += D.Asr * (GeraVap[D.ii][D.jj] * D.dt / (D.eps * D.Asr * rohl)) * D.eps * rohl;
					AgAcumu = AgAcumu - GeraVap[D.ii][D.jj] * D.dt;
                    NIVEL[D.ii] = NIVEL[D.ii] - (GeraVap[D.ii][D.jj] * D.dt / (D.eps * D.Asr * rohl));
					
					if (NIVEL[D.ii]>0) auxNIVEL = 1;
                }

            }
            ///////////////////////////////////////////////
            ////    INTERFACE NO ELEMENTO DE VOLUME    ////
            ///////////////////////////////////////////////
            else if( ((NIVEL[D.ii] / D.dL) >= aux) && ((NIVEL[D.ii]/D.dL) < (aux + 1)) && (NIVEL[D.ii] < D.AltLeito) && (D.mliq[D.ii][D.jj] > 0) )
            {
				
				CalculeiInterface[D.ii] = 1;
				
                if ((n == 0) && ((D.jj + nEsp) <= D.NumL))
                {
                    VazLiqEl = D.mliq[D.ii][D.jj];
                    for(int z = 0 ; z < nEsp ; z++)
                    {
                        D.mliq[D.ii][D.jj + z] = VazLiqEl * DistEspuma[z];
                    }
                }
                else
                {
                    D.Tflu[D.ii][D.jj] = D.Tflu[D.ii][D.jj-n];
                }

                if (D.ii == 0)
                {
                    Tvap = D.ToCoq;
                }
                else
                {
                    Tvap = D.Tflu[D.ii-1][D.jj+1];
                    if (Tvap < D.Tsat)
                    {
                        Tvap = D.Tsat;
                    }
                }
				
                HliqIn = PP.Entalpia(D.Tflu[D.ii][D.jj],D.P,0);
                rohl = PP.rohH2O(D.Tflu[D.ii][D.jj], 0);
				
                GeraVap[D.ii][D.jj] = 0;
                DT1 = D.Tsup[D.ii][D.jj] - D.Tflu[D.ii][D.jj];
                int contloop = 1;
                int contloop2 = 1;
                DVzOld = 0;
                do
                {					
                    L = Lvelho + (D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj])*D.dt / (D.eps * D.Asr * rohl);      //Aumento de nível devido a quantidade de água inserida [m]
                    D.Apt = 3 * (1 - D.eps) * D.pi * L * pow((Dmed / 2),2) / (D.Dp / 2);                        //Área superficial das partículas no elemento finito [m²]

                    Tsai = HLV.CalcTsai(D, PP, 0);

                    if (Tsai >= D.Tsat)
                    {   //houve ebulição
                        Tsai = D.Tsat;
                        DT2 = D.Tsup[D.ii][D.jj] - Tsai;
                        if ((((DT1 > 0) && (DT2 > 0)) || ((DT1 < 0) && (DT2 < 0))) && ((DT1 - DT2) != 0))
                        {
                            DTml = (DT1 - DT2) / log(DT1 / DT2);
                        }
                        else if ((DT1 - DT2) == 0)
                        {
                            DTml = D.Tsup[D.ii][D.jj] - D.Tsat;
                        }
                        else
                        {
                            DTml = 0;
                        }

                        Hsat = PP.Entalpia(D.Tsat, D.P, 0);
                        Qsensivel = D.mliq[D.ii][D.jj] * (Hsat - HliqIn);
                        DTe = D.Tsup[D.ii][D.jj] - D.Tsat;

                        if (DTe >= 120)
                        {   //ebulição em película
                            hebu = HLV.calchebu(D.Tsup[D.ii][D.jj], D, PP);
                            Qebu = hebu * D.Apt * DTe;
                            Qconv = Qsensivel;
                        }
                        else if (DTe >= 30)
                        {   //ebulição em regime de transição
                            hmin = HLV.calchebu((D.Tsat+120), D, PP);
                            a1 = (hmax - hmin) / (30 - 120);
                            b1 = hmax - a1 * 30;

                            hebu = (a1 * DTe + b1) * D.mhe;
                            Qebu = hebu * D.Apt * DTe;
                            Qconv = Qsensivel;
                        }
                        else if (DTe >= 5)
                        {   //ebulição nucleada
                            hmin = HLV.calch(D.Tflu[D.ii][D.jj], D.Tsat, (D.Tsat + 5), D.mliq[D.ii][D.jj], D, PP, 0);
                            a1 = (hmax - hmin) / (30 - 5);
                            b1 = hmax - a1 * 30;

                            hebu = (a1 * DTe + b1) * D.mhe;
                            Qebu = hebu * D.Apt * DTe;
                            Qconv = Qsensivel;
                        }
                        else
                        {   //ebulição livre
                            Qebu = 0;
                            Qconv = (HLV.h * D.Apt * DTml);
                        }

                        Qcoq = Qconv + Qebu;
                        hconv = Qcoq / (D.Apt * DTml);
                        GeraVap[D.ii][D.jj] = (Qcoq - Qsensivel) / PP.DHVapTsat;
                        if (GeraVap[D.ii][D.jj] > D.mliq[D.ii][D.jj])
                        {
                            GeraVap[D.ii][D.jj] = D.mliq[D.ii][D.jj];
                            Qebu = D.mliq[D.ii][D.jj] * PP.DHVapTsat;
                            Qcoq = Qconv + Qebu;
                            hconv = Qcoq / (D.Apt * DTml);
                        }
                    }
                    else
                    {   //não houve ebulição
                        hconv = HLV.h;
                        Qcoq = HLV.q;
                        GeraVap[D.ii][D.jj] = 0;
                    }

                    L2 = Lvelho + (D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj])*D.dt / (D.eps * D.Asr * rohl);
                    DIFF2 = abs((L2 - L) / L2);

                    //Evitando loop infinito
                    if (contloop >= D.nmax)
                    {
                        MsgResfriamento = "Atingiu numero maximo de iteracoes em RESFRIAMENTO";
                        if ((aviso2 == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgResfriamento << endl;
						aviso2 = 1;
						
                        if ((D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj]) == 0)
                        {
                            DIFF2 = 0;
                        }
                        else
                        {
                            if ((contloop2 >= 5) && (DVzOld > (D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj])))
                            {
                                DIFF2 = 0;
                            }
                            contloop2 ++;
                        }
                        DVzOld = (D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj]);
                    }
                    contloop ++;
                }while (DIFF2 > 0.00001);

                if (n == 0)
                {
                    NIVEL[D.ii+1] = NIVEL[D.ii] + (L - Lvelho);
					if (NIVEL[D.ii+1]>0) auxNIVEL = 1;
                }
                else
                {	
					auxNIVEL = NIVEL[D.ii+1];
                    NIVEL[D.ii+1] = NIVEL[D.ii+1] + (L - Lvelho);
					if (NIVEL[D.ii+1]>0) auxNIVEL = 1;
                }

                if ((n == nEsp) || ((D.jj + nEsp) > D.NumL))
                {
                    n = 0;
                }
                else
                {
                    n += 1;
                }

                rohv = PP.rohH2O(Tvap, 1);
                DeslocVol = (D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj]) * rohv / rohl;
                AgAcumu += D.Asr * (L - Lvelho) * D.eps * rohl;

                if (NIVEL[D.ii+1] > (D.jj * D.dL + D.dL))
                {
                    Lvelho = NIVEL[D.ii+1] - (D.jj * D.dL + D.dL);
                    D.Tflu[D.ii][D.jj+1] = Tsai;
                    //D.mliq[D.ii][D.jj+1] = D.mliq[D.ii][D.jj] - GeraVap[D.ii][D.jj];
                }
                else
                {
                    Lvelho = L;
                    D.Tflu[D.ii][D.jj+1] = mistura_correntes(D.Tsat, Tvap, GeraVap[D.ii][D.jj], DeslocVol, D, PP, 1);
                }
                D.mvap[D.ii][D.jj+1] = D.mvap[D.ii][D.jj] + GeraVap[D.ii][D.jj] + DeslocVol;

                HLV.h = hconv;
                HLV.Tbulk = (D.Tflu[D.ii][D.jj] + D.Tflu[D.ii][D.jj+1])/2;
                HLV.q = Qcoq;

                D.Tsup[D.ii+1][D.jj] = KLV.CalcTsup(D, PP, HLV);
                D.Tcoq[D.ii+1][D.jj] = KLV.CalcTcoq(D, PP, HLV);
				
				/*if (D.ii == 4201 && D.jj == 59)
				{
					cout << endl << "deu ruim - 2" << endl;
					cout << "D.ii = " << D.ii << " D.jj = " << D.jj << endl;
					cout << "GeraVap[D.ii][D.jj] = " << GeraVap[D.ii][D.jj] << endl;
					cout << "D.mliq[D.ii][D.jj+1] = " << D.mliq[D.ii][D.jj+1] << endl;
					cout << "D.mvap[D.ii][D.jj+1] = " << D.mvap[D.ii][D.jj+1] << endl;
					cout << "D.Tflu[D.ii][D.jj+1] = " << D.Tflu[D.ii][D.jj+1] << endl;
					cout << "D.Tsup[D.ii+1][D.jj] = " << D.Tsup[D.ii+1][D.jj] << endl;
					cout << "D.Tcoq[D.ii+1][D.jj] = " << D.Tcoq[D.ii+1][D.jj] << endl;
					
					cout << "hconv = " << hconv << endl;
					cout << "HLV.Tbulk = " << HLV.Tbulk << endl;
					cout << "HLV.q = " << HLV.q << endl;
					
					system("PAUSE");
				}*/
				
            }
            /////////////////////////////////////////
            ////    ELEMENTO NÃO RECEBEU ÁGUA    ////
            /////////////////////////////////////////
            else
            {	
                GeraVap[D.ii][D.jj] = 0;
                D.mliq[D.ii][D.jj+1] = 0;
                D.mvap[D.ii][D.jj+1] = D.mvap[D.ii][D.jj];

                D.Tflu[D.ii][D.jj+1] = HLV.CalcTsai(D, PP, 1);
				
                D.Tsup[D.ii+1][D.jj] = KLV.CalcTsup(D, PP, HLV);
                D.Tcoq[D.ii+1][D.jj] = KLV.CalcTcoq(D, PP, HLV);
				
				/*if (D.ii == 4201 && D.jj == 59)
				{
					cout << "deu ruim - 3" << endl;
					cout << "D.ii = " << D.ii << " D.jj = " << D.jj << endl;
					cout << "GeraVap[D.ii][D.jj] = " << GeraVap[D.ii][D.jj] << endl;
					cout << "D.mliq[D.ii][D.jj+1] = " << D.mliq[D.ii][D.jj+1] << endl;
					cout << "D.mvap[D.ii][D.jj+1] = " << D.mvap[D.ii][D.jj+1] << endl;
					cout << "D.Tflu[D.ii][D.jj+1] = " << D.Tflu[D.ii][D.jj+1] << endl;
					cout << "D.Tsup[D.ii+1][D.jj] = " << D.Tsup[D.ii+1][D.jj] << endl;
					cout << "D.Tcoq[D.ii+1][D.jj] = " << D.Tcoq[D.ii+1][D.jj] << endl;
					
					system("PAUSE");
				}*/
				
            }

            VaporGeradoTotal += GeraVap[D.ii][D.jj] * D.dt;
            if (D.mvap[D.ii][D.jj+1] > PicoVazao)
            {
                PicoVazao = D.mvap[D.ii][D.jj+1];
            }
			
			/*if ((D.ii == 4200) )
			{
				cout << endl << "D.ii = " << D.ii << " D.jj = " << D.jj << endl;
				if(D.jj == 0) cout << "D.Tflu[" << D.ii << "][" << D.jj << "] = " <<  D.Tflu[D.ii][D.jj] << endl;
				cout << "D.Tflu[" << D.ii << "][" << D.jj+1 << "] = " <<  D.Tflu[D.ii][D.jj+1] << endl;
				cout << "D.Tsup[" << D.ii+1 << "][" << D.jj << "] = " <<  D.Tsup[D.ii+1][D.jj] << endl;
				cout << "D.Tcoq[" << D.ii+1 << "][" << D.jj << "] = " <<  D.Tcoq[D.ii+1][D.jj] << endl;
				system("PAUSE");
			}*/
			
        }//Next jj
		
		if ((NIVEL[D.ii+1] == 0) && (NIVEL[D.ii] > 0) && (auxNIVEL == 1))
		{	
			NIVEL[D.ii+1] = NIVEL[D.ii];
		}			

        if (NIVEL[D.ii] >= D.AltLeito)
        {   //Nivel ultrapaddou a altura do coque no leito
            rohl = PP.rohH2O(D.Tsat, 0);
            L = D.mliq[D.ii][D.NumL-1] * D.dt / (D.Asr * rohl);
            NIVEL[D.ii+1] = NIVEL[D.ii] + L;
            D.Tflu[D.ii][D.jj] = D.Tsat;
            rohv = PP.rohH2O(D.Tsat, 1);
            D.mvap[D.ii][D.NumL] += D.mliq[D.ii][D.NumL] * rohv/rohl;
            AgAcumuAcimaLeito += D.Asr * L * rohl;
			
			if (NIVEL[D.ii+1] > 0) auxNIVEL = 1;			
        }
		
		//cout << "D.Tflu[" << D.ii << "][" << D.jj << "] = " <<  D.Tflu[D.ii][D.jj] << endl;
		
    }//Next ii
	
    AgAcumu += AgAcumuAcimaLeito;

	//system("PAUSE");
}

void RESFRIAMENTO::espuma(DADOS &D)
{
    nEsp = D.Le/D.dL;
    double aux1 = D.Le/D.dL - nEsp;
    if (aux1 > 0.0)
    {
        nEsp += 1;
    }

    DistEspuma.resize(nEsp);

    double fla = 0;
    for (int i = 1; i <= nEsp ; i++)
    {
        fla += 1 / (pow(i,D.zd));
    }
    double AA = 1/fla;

    for (int i = 1; i <= nEsp ; i++)
    {
        DistEspuma[i-1] = AA/(pow(i,D.zd));
    }
}

double RESFRIAMENTO::mistura_correntes(double T1, double T2, double m1, double m2, DADOS &D, PROPRIEDADES &PP, int estado)
{
    double H1;
    double H2;
    double Hmx;
    double Hmx2;
    double Hmx2dT;
    double Tmx;
    double fT;
    double fTdT;
    double flT;
    double dTempera = 0.001;

    if ((estado == 1) && (T1 < D.Tsat))
    {
        T1 = D.Tsat;
    }
    if ((estado == 1) && (T2 < D.Tsat))
    {
        T2 = D.Tsat;
    }

    H1 = PP.Entalpia(T1,D.P,estado);
    H2 = PP.Entalpia(T2,D.P,estado);
    Hmx = (m1 * H1 + m2 * H2)/(m1 + m2);

    //Estimativa inicial de temperatura
    Tmx = (m1 * T1 + m2 * T2)/(m1 + m2);
    if ((estado == 1) && (Tmx < D.Tsat))
    {
        Tmx = D.Tsat;
    }

    int contloop = 1;
    do
    {
        Hmx2 = PP.Entalpia(Tmx,D.P,estado);
        fT = Hmx - Hmx2;
        Hmx2dT = PP.Entalpia((Tmx + dTempera),D.P,estado);
        fTdT = Hmx - Hmx2dT;
        flT = (fTdT - fT) / dTempera;
        Tmx = Tmx - (fT / flT);

        //evitando loop infinito
        if (contloop == D.nmax)
        {
            MsgResfriamento = "Atingiu numero maximo de iteraçlões em mistura_correntes";
            if ((aviso == 0) && (D.PSO == 0)) cout << D.ii << " " << D.jj << " " << MsgResfriamento << endl;
			aviso = 1;
            fT = 0;
        }
        contloop += 1;
    }while (fT > 0.00000001);

    return Tmx;

    }



//RESFRIAMENTO::~RESFRIAMENTO()
//{
    //dtor
//}
