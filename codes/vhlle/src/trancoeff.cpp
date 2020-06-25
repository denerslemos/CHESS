#include "trancoeff.h"
#include "eos.h"
#include "inc.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, EoS *_eos)
{
  etaS = _etaS;
  zetaS = _zetaS;
  eos = _eos;
}

void TransportCoeff::getEta(double e, double T, double &_etaS, double &_zetaS)
{
/*  double Tc = 184.;
  double TT = T/Tc;
  double aa=0.15;
  double bb=0.14;
  double cc=0.66;
  double deltaa=5.1;
  double betaa = 33/(12*3.14159);
  double gammaa = 1.6;
  double alphaa = (1/betaa)*((pow(cc*TT,2)-1)/(pow(cc*TT,2) * log(pow(cc*TT,2))));
  double zeta = aa/pow(alphaa,gammaa) + bb/pow(TT,deltaa);
  _etaS = etaS*zeta;
*/
  _etaS = etaS;
  _zetaS = zetaS * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.);
}


void TransportCoeff::getTau(double T, double &_taupi, double &_tauPi)
{
  if (T > 0.)
    _taupi = std::max(5. / 5.068 * etaS / T, 0.003);
  else
    _taupi = 0.1;
  if (T > 0.)
    _tauPi = std::max(5. / 5.068 * (1. / 4. / C_PI) / T, 0.005);
  else
    _tauPi = 0.1;
}
