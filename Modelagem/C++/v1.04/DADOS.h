#ifndef DADOS_H
#define DADOS_H

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class DADOS
{
    public:

    //Utilidade pública
        long double pi;                     //numero pi
        double Asr;                         //Área da secção reta do leito [m²]
        double Apt;                         //Área das partículas no elemento finito [m²]
        double Tsat;                        //temperatura de saturação na presão do resfriamento [K]
        double Mcoq;                        //Massa de coque no elemento [kg]
        int ii;                             //contador de tempo
        int jj;                             //contador de espaço
        int NumL;                           //numero de elementos de espaço
        int NumT;                           //numero de elementos no tempo
		int seed;							//semente do numero randomico

        vector<vector<double> > Tflu;        //Temperatura do fluido na entrada do elemento finito num dado instante do tempo [K]
        vector<vector<double> > Tcoq;        //Temperatura média do interior do coque no elemento finito num dado intante do tempo [K]
        vector<vector<double> > Tsup;        //Temperatura média da superfície do coque no elemento finito num dado intante do tempo [K]
        vector<vector<double> > mliq;        //vazão mássica de liquido na entrada do elemento de volume [kg/s]
        vector<vector<double> > mvap;        //vazão mássica de vapor na entrada do elemento de volume [kg/s]

    //DADOS DE GEOMETRIA
        double Dc;                          //diâmetro do leito[m]
        double AltLeito;                    //Altura do Leito de coque [m]
        double Dcone;                       //dimâmetro menor do cone [m]
        double hcone;                       //altura do cone [m]

    //DADOS DE PROCESSO
        double P;                           //Pressão do leito [bar abs]
        double ToCoq;                       //Temperatura inicial do coque [K]
        double Ttot;                        //Duração da purga com vapor [s]

        vector<vector<double> > CurvResf;    //matriz de curva de vazão mássica de água ou vapor para o resfriamento ou purga. linhas: patamares - colunas: vazãos [kg/h], tempo acumulado[min], temperatura [°C]

    //PROPRIEDADES DO COQUE
        double rohcoq;                      //massa específica do coque [kg/m³]
        double cpcoq;                       //capacidade calorífica do coque [kJ/kg.K]
        double Kcoq;                        //condutividade térmica do coque [kW/m.K]
        double eps;                         //porosidade do leito
        double emiss;                       //emissividade do coque

    //PARÂMETROS DE AJUSTE DO MODELO
        double psi;                         //tortuosidade do leito
        double Dp;                          //Diâmetro médio das particulas de coque no leito [m]
        double mhl;                         //fator de ajuste do coeficiente de convecção livre
        double mhc;                         //fator de ajuste do coeficiente de convecção forçada
        double mhe;                         //fator de ajuste do coeficiente de convecção ebulição pelicula
        double Le;                          //comprimento da interface/espuma [m]
        double zd;                          //fator para inclinação da reta

    //PARâMETROS NUMÉRICOS DO MODELO
        double dL;                          //Altura do elemento finito de volume [m]
        double dt;                          //Passo de integração no tempo [s]
        int NamsT;                          //numero de amostras de tempo para exibição do resultado
        int nmax;                           //numero máximo de iterações

    //PARA O ENXAME DE PARTÍCULAS (PSO)
        bool PSO;                           //escolha do usuário para realização do PSO
        int Iteracoes;                      //numero de iterações do PSO
        int Particulas;                     //numero de particulas para o PSO
        int Parametros;                     //numero de parâmetros para o PSO
		int NddsCrvExp;						//número de dados da curva experimental de temperatura na solda
		double AgVapExp;					//Total de água vaporizada experimental (medido no blowdown)

        vector<double> xMax;                //valores(posições) de máximo dos parÂmetros
        vector<double> xMin;                //valores(posições) de mínimo dos parâmetros
        vector<vector<double> > ResfExp;	//dados experimentais de temperatura [°C] vs. tempo [min]
		vector<int> auxPar;				//vetor auxiliar para os parametros

        DADOS();

};

#endif // DADOS_H
