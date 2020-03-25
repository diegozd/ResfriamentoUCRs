#ifndef DADOS_H
#define DADOS_H

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class DADOS
{
    public:

    //Utilidade p�blica
        long double pi;                     //numero pi
        double Asr;                         //�rea da sec��o reta do leito [m�]
        double Apt;                         //�rea das part�culas no elemento finito [m�]
        double Tsat;                        //temperatura de satura��o na pres�o do resfriamento [K]
        double Mcoq;                        //Massa de coque no elemento [kg]
        int ii;                             //contador de tempo
        int jj;                             //contador de espa�o
        int NumL;                           //numero de elementos de espa�o
        int NumT;                           //numero de elementos no tempo
		int seed;							//semente do numero randomico

        vector<vector<double> > Tflu;        //Temperatura do fluido na entrada do elemento finito num dado instante do tempo [K]
        vector<vector<double> > Tcoq;        //Temperatura m�dia do interior do coque no elemento finito num dado intante do tempo [K]
        vector<vector<double> > Tsup;        //Temperatura m�dia da superf�cie do coque no elemento finito num dado intante do tempo [K]
        vector<vector<double> > mliq;        //vaz�o m�ssica de liquido na entrada do elemento de volume [kg/s]
        vector<vector<double> > mvap;        //vaz�o m�ssica de vapor na entrada do elemento de volume [kg/s]

    //DADOS DE GEOMETRIA
        double Dc;                          //di�metro do leito[m]
        double AltLeito;                    //Altura do Leito de coque [m]
        double Dcone;                       //dim�metro menor do cone [m]
        double hcone;                       //altura do cone [m]

    //DADOS DE PROCESSO
        double P;                           //Press�o do leito [bar abs]
        double ToCoq;                       //Temperatura inicial do coque [K]
        double Ttot;                        //Dura��o da purga com vapor [s]

        vector<vector<double> > CurvResf;    //matriz de curva de vaz�o m�ssica de �gua ou vapor para o resfriamento ou purga. linhas: patamares - colunas: vaz�os [kg/h], tempo acumulado[min], temperatura [�C]

    //PROPRIEDADES DO COQUE
        double rohcoq;                      //massa espec�fica do coque [kg/m�]
        double cpcoq;                       //capacidade calor�fica do coque [kJ/kg.K]
        double Kcoq;                        //condutividade t�rmica do coque [kW/m.K]
        double eps;                         //porosidade do leito
        double emiss;                       //emissividade do coque

    //PAR�METROS DE AJUSTE DO MODELO
        double psi;                         //tortuosidade do leito
        double Dp;                          //Di�metro m�dio das particulas de coque no leito [m]
        double mhl;                         //fator de ajuste do coeficiente de convec��o livre
        double mhc;                         //fator de ajuste do coeficiente de convec��o for�ada
        double mhe;                         //fator de ajuste do coeficiente de convec��o ebuli��o pelicula
        double Le;                          //comprimento da interface/espuma [m]
        double zd;                          //fator para inclina��o da reta

    //PAR�METROS NUM�RICOS DO MODELO
        double dL;                          //Altura do elemento finito de volume [m]
        double dt;                          //Passo de integra��o no tempo [s]
        int NamsT;                          //numero de amostras de tempo para exibi��o do resultado
        int nmax;                           //numero m�ximo de itera��es

    //PARA O ENXAME DE PART�CULAS (PSO)
        bool PSO;                           //escolha do usu�rio para realiza��o do PSO
        int Iteracoes;                      //numero de itera��es do PSO
        int Particulas;                     //numero de particulas para o PSO
        int Parametros;                     //numero de par�metros para o PSO
		int NddsCrvExp;						//n�mero de dados da curva experimental de temperatura na solda
		double AgVapExp;					//Total de �gua vaporizada experimental (medido no blowdown)

        vector<double> xMax;                //valores(posi��es) de m�ximo dos par�metros
        vector<double> xMin;                //valores(posi��es) de m�nimo dos par�metros
        vector<vector<double> > ResfExp;	//dados experimentais de temperatura [�C] vs. tempo [min]
		vector<int> auxPar;				//vetor auxiliar para os parametros

        DADOS();

};

#endif // DADOS_H
