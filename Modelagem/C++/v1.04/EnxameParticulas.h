#ifndef EnxameParticulas_H
#define EnxameParticulas_H

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>

#include "DADOS.h"
#include "PROPRIEDADES.h"
#include "RESFRIAMENTO.h"
#include "Random.h"

using namespace std;

class EnxameParticulas
{
    public:

        double c1;
        double c2;
        double wIni;
        double wFim;
        double w;
        double vPg ;
        double vPgN;

        vector<double> xMin;
        vector<double> xMax;
        vector<double> Pg;
        vector<double> FobjValor;
        vector<double> vPl;
        vector<double> vMax;
        vector<vector<double> > X;
        vector<vector<double> > Pl;
        vector<vector<double> > V;

        vector<vector<double> > MatrizCompara;

        EnxameParticulas(DADOS&);
        double FuncaoObjetivo(vector<double>);
        void AvaliarFobj(DADOS&);
        void AvaliarSolucao(DADOS&);
        void AtualizarPosicao(DADOS&);
        void Otimizar(DADOS&);


    private:
        int NumDados;
        int solda;
        double aux;


};
#endif
