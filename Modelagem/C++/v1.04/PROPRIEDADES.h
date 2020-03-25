#ifndef PROPRIEDADES_H
#define PROPRIEDADES_H

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#include "DADOS.h"

using namespace std;

class PROPRIEDADES
{
    public:

        double Tsat;        //temperatura de satura��o na pres�o do resfriamento [K]
        double DHVapTsat;   //Entalpia de vaporiza��o na temperatura de satura��o [kJ/kg]

        string MsgPropriedades;

        PROPRIEDADES(DADOS&);

        double rohH2O(double, int);
        double cpH2O(double, int);
        double KH2O(double, int);
        double muH2O(double, int);
        double sigH2O(double);
        double betaH2O(double, int);
        double DHvapH2O(double);
        double Entalpia(double, double, int);

    private:
        double aux;
        int vai;

};
#endif // PROPRIEDADES_H
