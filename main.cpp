#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
class polinom
{
    int grad; //gradul polinomului
    double* coef; //coeficientii
public:
    polinom();
    polinom(int g); //polinom de grad g, cu coeficienti initial nuli
    polinom(const polinom &p);//constructor de copiere
    ~polinom();
    double calcul_valoare(double x);
    bool operator==(const polinom& p);
    friend istream& operator>>(istream& i, polinom& p);
    friend ostream& operator<<(ostream& o, const polinom& p);
    friend polinom operator+(const polinom& p1,const polinom& p2); //nu e metoda a clasei, deoarece numarul de coef difera
    friend polinom operator-(const polinom& p1,const polinom& p2);
    friend polinom operator*(const polinom& p1,const polinom& p2);
    friend class pereche;
};

class pereche
{
    double r;
    polinom p;
public:
    pereche();
    bool verificare_radacina();
    bool operator==(const pereche& p);
    friend istream& operator>>(istream& i, pereche& pr);
    friend ostream& operator<<(ostream& i, const pereche& pr);
};

polinom::polinom()
{
    grad = -1; //conventie
    coef = new double;
    *coef = 0;
}

polinom::polinom(int g)
{
    grad = g;
    if(grad == -1) //polinom nul, apare la diferenta, cand doua polinoame sunt identice
    {
        coef = new double;
        *coef = 0;
    }
    else
    if(grad == 0) //cand gradul este 0, deci polinomul este format doar din termenul liber
    {
        coef = new double;
        *coef = 0; //initial 0, pentru a fi completat
    }
    else //grad >= 1
    {
        coef = new double[grad+1];
        for(int i = 0; i <= grad; i++)
        coef[i] = 0;
    }
}

polinom::polinom(const polinom &p)
{
    grad = p.grad;
    coef = new double[grad+1];
    for(int i = 0; i <= grad; i++)
        coef[i] = p.coef[i];
}

polinom::~polinom()
{
    delete []coef; //eliberez zona de memorie asociata vectorului de coeficienti
}

double polinom::calcul_valoare(double x)
{
    if(grad == -1) //polinomul nul
        return 0;
    if(grad == 0) //polinomul format din termenul liber
        return *coef;
    double S = 0;
    for(int i = 0; i <= grad; i++)
        S+= pow(x,i)*coef[i];
    return S;
}

bool polinom::operator==(const polinom& p)
 {
     if(grad != p.grad) //daca gradele sunt diferite, polinoamele sunt diferite
        return 0;
     for(int i = 0; i <= grad; i++)
        if(coef[i] != p.coef[i]) //compar coeficient cu coeficient
        return 0;
     return 1;
 }

istream& operator>>(istream& i,polinom& p)
{
    i>>p.grad;
    delete []p.coef; //sterg zona de memorie alocata in constructor
    p.coef = new double[p.grad+1]; //aloc o noua zona de memorie conform gradului
    for(int k = 0; k <= p.grad; k++)
        i>>p.coef[k];
    return i;
}

ostream& operator<<(ostream& o,const polinom& p)
{
    if(p.grad == -1)
        o<<0<<" ";
    else
    {
        polinom nul = polinom(p.grad); //polinomul cu coef nuli de grad egal cu cel de afisat
        if(nul == p) // daca polinomul are toti coeficientii nuli
        {
            o<<0<<" ";
            return o;
        }
        int start = 0;
        while(p.coef[start] == 0)
            start++; //pana la k, avem coeficienti nuli, deci nu afisam nimic
        if(start == 0) //daca termenul liber este nenul
        o<<p.coef[start]<<" ";
        else
        o<<p.coef[start]<<"x^"<<start<<" ";// primul monom nenul
    for(int k = start+1 ; k <= p.grad; k++)
        if(p.coef[k] == 1) //daca un coeficient este 1, il omit
        o<<"+"<<"x^"<<k<<" ";
        else
        if(p.coef[k] == -1) //daca un coeficient este -1, las doar semnul -
        o<<"-x^"<<k<<" ";
        else
        if(p.coef[k] > 0) //daca un coeficient este pozitiv, pun semnul + inainte
        o<<" +"<<p.coef[k]<<"x^"<<k<<" ";
        else
        if(p.coef[k] < 0) //daca un coeficient este negativ, semnul - apare oricum
        o<<p.coef[k]<<"x^"<<k<<" ";
    }
    return o;
}

polinom operator+(const polinom& p1,const polinom& p2)
{
    int grad_S; //gradul polinomului suma este maximul dintre gradele celuilalt
    if(p1.grad > p2.grad)
        grad_S = p1.grad;
    else
        grad_S = p2.grad;
    polinom S(grad_S); //creez un polinom de grad maxim cu coeficienti initial nuli
    int i = 0;
    while(i <= p1.grad && i <= p2.grad) //adun coeficient cu coeficient pana ajung la finalul polinomului de grad minim
    {
        S.coef[i] = p1.coef[i] + p2.coef[i];
        i++;
    }
    while(i <= p1.grad) //daca primul polinom e de grad mai mare decat al doilea
    {
        S.coef[i] = p1.coef[i];
        i++;
    }
     while(i <= p2.grad) //daca al doilea polinom e de grad mai mare decat primul
    {
        S.coef[i] = p2.coef[i];
        i++;
    }
    return S;
}

polinom operator-(const polinom& p1,const polinom& p2)
{
    int grad_D; //gradul polinomului diferenta este initial maximul dintre gradele celuilalt
    if(p1.grad == p2.grad)
    {
        grad_D = p1.grad;
        while(p1.coef[grad_D] == p2.coef[grad_D] && grad_D >=0) //daca se reduc termenii de grad maxim, scade gradul dif
        grad_D--;
        if(grad_D == -1) //daca gradul este -1, cele doua polinoame sunt identice, deci diferenta este polinomul nul
            {
                polinom D(-1);
                return D;
            }
    }
    else //daca cele doua polinoame nu au acelasi grad, gradul dif este  maximul dintre gradele celuilalt
    {
        if(p1.grad > p2.grad)
        grad_D = p1.grad;
    else
        grad_D = p2.grad;
    }
    polinom D(grad_D); //creez polinomul de grad maxim cu coeficienti initial nuli
    int i = 0;
    while(i <= p1.grad && i <= p2.grad && i <= grad_D) //scad coeficient cu coeficient pana ajung la finalul polinomului de grad minim, sau daca sunt egale, pana la gradul la care se reduc (grad_D)
    {
        D.coef[i] = p1.coef[i] - p2.coef[i];
        i++;
    }
    while(i <= p1.grad && i <= grad_D) //daca primul polinom e de grad mai mare decat al doilea
    {
        D.coef[i] = p1.coef[i];
        i++;
    }
     while(i <= p2.grad && i <= grad_D) //daca al doilea polinom e de grad mai mare decat al doilea
    {
        D.coef[i] = -p2.coef[i];
        i++;
    }
    return D;
}

polinom operator*(const polinom& p1,const polinom& p2)
{
    int grad_P = p1.grad + p2.grad; //gradul polinomului produs este suma gradelor celor doua polinoame
    polinom P(grad_P); //creez un polinom de grad corespunzator cu coeficienti initial nuli
    for(int i = 0; i <= p1.grad; i++)
        for(int j = 0; j <= p2.grad; j++)
        P.coef[i+j] += p1.coef[i]*p2.coef[j]; //adun pe pozitia corespunzatoare produsul fiecarui coeficient din primul polinom cu fiecare coeficient din al doilea
    return P;
}

pereche::pereche()
{
    r = 0;
    polinom p(-1);
}

bool pereche::operator==(const pereche& pr)
{
  if(r != pr.r)
    return false;
  if(p == pr.r)
    return true;
  return false;
}

bool pereche::verificare_radacina()
{
    if(p.calcul_valoare(r) == 0) //daca polinomul p ia valoarea 0 pentru valoarea reala r, atunci e radacina
        return true;
    return false;
}

istream& operator>>(istream& i, pereche& pr)
{
    i>>pr.p;//citesc polinomul
    i>>pr.r;//citesc numarul real
    return i;
}

ostream& operator<<(ostream& i,const pereche& pr)
{
    i<<"polinomul:"<<pr.p<<"\t";
    i<<"numarul real:"<<pr.r;
    return i;
}


int main()
{
    //citire din fisier si scriere in fisier
    ifstream f;
    f.open("input.txt");
    ofstream g;
    g.open("output.txt");
    int n;
    f>>n;
    polinom vp[n];
    for(int i = 0; i < n; i++)
        f>>vp[i];
    for(int i = 0; i < n; i++) //afisez
        g<<"Polinomul "<<i+1<<": "<<vp[i]<<endl;
    g<<endl;
    for(int i = 0; i < n; i++)
        g<<"Suma dintre polinomul 1 si polinomul "<<i+1<<": "<<vp[0]+vp[i]<<endl;
    g<<endl;
    for(int i = 0; i < n; i++)
        g<<"Diferenta dintre polinomul 1 si polinomul "<<i+1<<": "<<vp[0]-vp[i]<<"."<<endl;
    g<<endl;
    for(int i = 0; i < n; i++)
        g<<"Produsul dintre polinomul 1 si polinomul "<<i+1<<": "<<vp[0]*vp[i]<<"."<<endl;
    g<<endl;

    int n_p;
    f>>n_p;
    pereche vpr[n_p];
    for(int i = 0; i < n_p; i++)
        f>>vpr[i];
    for(int i = 0; i < n_p; i++) //afisez perechile
    {
        g<<"Perechea "<<i+1<<": "<<vpr[i]<<endl;
        if(vpr[i].verificare_radacina() == true)
            g<<"Numarul real este radacina pentru polinom."<<endl;
        else
            g<<"Numarul real nu este radacina pentru polinom."<<endl;
    }
    for(int i = 0; i < n_p; i++)
      if(vpr[i] == vpr[0])
      g<<"Perechea 1 este egala cu perechea "<<i+1<<". "<<endl;
    else
      g<<"Perechea 1 nu este egala cu perechea "<<i+1<<". "<<endl;
    f.close();
    g.close();

  /*
    // citire de la tastatura si scriere pe ecran
    polinom p;
    cout<<"Introduceti polinomul (grad si coeficientii, pornind de la termenul liber pana la coeficientul de grad maxim):"<<endl;
    cin>>p;
    cout<<"Polinomul la puterea a treia este: "<<p*p*p<<endl;
*/
    return 0;
}
