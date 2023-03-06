#include <iostream>
#include <string>
#include <ctime>
#include <string>
#include <vector>


using namespace std;

struct atom_param{
    string type;
    double radius;
    double epsilon;
    double mass;
};

int main(int argc,char **argv){
    atom_param ap;
    vector<atom_param> gaff_parameters;

    ap.type = hc;
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = ha;
    ap.radius = 1.4735;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = hn;
    ap.radius = 0.6210;
    ap.epsilon = 0.0100;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = ho;
    ap.radius = 0.3019;
    ap.epsilon = 0.0047;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = hs;
    ap.radius = 0.6112;
    ap.epsilon = 0.0124;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = hp;
    ap.radius = 0.6031;
    ap.epsilon = 0.0144;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = o;
    ap.radius = 1.7107;
    ap.epsilon = 0.1463;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = os;
    ap.radius = 1.7713;
    ap.epsilon = 0.0726;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = oh;
    ap.radius = 1.8200;
    ap.epsilon = 0.0930;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = ow;
    ap.radius = 1.7683;
    ap.epsilon = 0.1520;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = c3;
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = c2;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = c1;
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = n;
    ap.radius = 1.7852;
    ap.epsilon = 0.1636;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = s;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = p2;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = f;
    ap.radius = 1.7029;
    ap.epsilon = 0.0832;
    ap.mass = 19.00;
    gaff_parameters.append(ap);

    ap.type = cl;
    ap.radius = 1.9452;
    ap.epsilon = 0.2638;
    ap.mass = 35.45;
    gaff_parameters.append(ap);

    ap.type = br;
    ap.radius = 2.0275;
    ap.epsilon = 0.3932;
    ap.mass = 79.90;
    gaff_parameters.append(ap);

    ap.type = i;
    ap.radius = 2.1558;
    ap.epsilon = 0.4955;
    ap.mass = 126.9;
    gaff_parameters.append(ap);

    ap.type = n1;
    ap.radius = 1.8372;
    ap.epsilon = 0.1098;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n2;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n3;
    ap.radius = 1.8886;
    ap.epsilon = 0.0858;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = na;
    ap.radius = 1.7992;
    ap.epsilon = 0.2042;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nh;
    ap.radius = 1.7903;
    ap.epsilon = 0.2150;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n+;
    ap.radius = 1.6028;
    ap.epsilon = 0.7828;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n9;
    ap.radius = 2.2700;
    ap.epsilon = 0.0095;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = h1;
    ap.radius = 1.3593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = h2;
    ap.radius = 1.2593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = h3;
    ap.radius = 1.1593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = hx;
    ap.radius = 1.0593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = h4;
    ap.radius = 1.4235;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = h5;
    ap.radius = 1.3735;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = cx;
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cy;
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = c;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cs;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = ca;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cc;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cd;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = ce;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cf;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cp;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cq;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cz;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cg;
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = ch;
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cu;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = cv;
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = nb;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nc;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nd;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = ne;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nf;
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = no;
    ap.radius = 1.8886;
    ap.epsilon = 0.0858;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n7;
    ap.radius = 1.9686;
    ap.epsilon = 0.0522;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n8;
    ap.radius = 2.0486;
    ap.epsilon = 0.0323;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = n4;
    ap.radius = 1.4028;
    ap.epsilon = 3.8748;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nx;
    ap.radius = 1.4528;
    ap.epsilon = 2.5453;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = ny;
    ap.radius = 1.5028;
    ap.epsilon = 1.6959;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nz;
    ap.radius = 1.5528;
    ap.epsilon = 1.1450;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = ns;
    ap.radius = 1.8352;
    ap.epsilon = 0.1174;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nt;
    ap.radius = 1.8852;
    ap.epsilon = 0.0851;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nu;
    ap.radius = 1.8403;
    ap.epsilon = 0.1545;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = nv;
    ap.radius = 1.8903;
    ap.epsilon = 0.1120;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = s2;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = s4;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = s6;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = sx;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = sy;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = sh;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = ss;
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = p3;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = p4;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = p5;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = pb;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = px;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = py;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = pc;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = pd;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = pe;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = pf;
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = H;
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HO;
    ap.radius = 0.0000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HS;
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HC;
    ap.radius = 1.4870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = H1;
    ap.radius = 1.3870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = H2;
    ap.radius = 1.2870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = H3;
    ap.radius = 1.1870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HP;
    ap.radius = 1.1000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HA;
    ap.radius = 1.4590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = H4;
    ap.radius = 1.4090;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = H5;
    ap.radius = 1.3590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HW;
    ap.radius = 0.0000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = HZ;
    ap.radius = 1.4590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = O;
    ap.radius = 1.6612;
    ap.epsilon = 0.2100;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = O2;
    ap.radius = 1.6612;
    ap.epsilon = 0.2100;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = OW;
    ap.radius = 1.7683;
    ap.epsilon = 0.1520;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = OH;
    ap.radius = 1.7210;
    ap.epsilon = 0.2104;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = OS;
    ap.radius = 1.6837;
    ap.epsilon = 0.1700;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = OP;
    ap.radius = 1.8500;
    ap.epsilon = 0.1700;
    ap.mass = 16.00;
    gaff_parameters.append(ap);

    ap.type = C*;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CI;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = C5;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = C4;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CT;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CX;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = C;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = N;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = N3;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = S;
    ap.radius = 2.0000;
    ap.epsilon = 0.2500;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = SH;
    ap.radius = 2.0000;
    ap.epsilon = 0.2500;
    ap.mass = 32.06;
    gaff_parameters.append(ap);

    ap.type = P;
    ap.radius = 2.1000;
    ap.epsilon = 0.2000;
    ap.mass = 30.97;
    gaff_parameters.append(ap);

    ap.type = MG;
    ap.radius = 0.7926;
    ap.epsilon = 0.8947;
    ap.mass = 24.305;
    gaff_parameters.append(ap);

    ap.type = C0;
    ap.radius = 1.7131;
    ap.epsilon = 0.45979;
    ap.mass = 40.08;
    gaff_parameters.append(ap);

    ap.type = Zn;
    ap.radius = 1.10;
    ap.epsilon = 0.0125;
    ap.mass = 65.4;
    gaff_parameters.append(ap);

    ap.type = F;
    ap.radius = 1.75;
    ap.epsilon = 0.061;
    ap.mass = 19.00;
    gaff_parameters.append(ap);

    ap.type = Cl;
    ap.radius = 1.948;
    ap.epsilon = 0.265;
    ap.mass = 35.45;
    gaff_parameters.append(ap);

    ap.type = Br;
    ap.radius = 2.22;
    ap.epsilon = 0.320;
    ap.mass = 79.90;
    gaff_parameters.append(ap);

    ap.type = I;
    ap.radius = 2.35;
    ap.epsilon = 0.40;
    ap.mass = 126.9;
    gaff_parameters.append(ap);

    ap.type = EP;
    ap.radius = 0.00;
    ap.epsilon = 0.0000;
    ap.mass = 0.00;
    gaff_parameters.append(ap);

    ap.type = 2C;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = 3C;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = C8;
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CO;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CA;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CB;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CC;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CD;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CK;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CM;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CQ;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CV;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CW;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CR;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CN;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CY;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CZ;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CP;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = CS;
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.append(ap);

    ap.type = N2;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = NA;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = NB;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = NC;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = NT;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = NY;
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.append(ap);

    ap.type = Fe;
    ap.radius = 1.200;
    ap.epsilon = 0.0500;
    ap.mass = 55.00;
    gaff_parameters.append(ap);

    ap.type = Cu;
    ap.radius = 2.200;
    ap.epsilon = 0.2000;
    ap.mass = 63.55;
    gaff_parameters.append(ap);

    ap.type = Ca;
    ap.radius = 1.790;
    ap.epsilon = 0.0140;
    ap.mass = 40.00;
    gaff_parameters.append(ap);

    ap.type = Si;
    ap.radius = 2.220;
    ap.epsilon = 0.3200;
    ap.mass = 28.00;
    gaff_parameters.append(ap);

    ap.type = hw;
    ap.radius = 0.000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.append(ap);

    ap.type = K;
    ap.radius = 2.658;
    ap.epsilon = 0.00033;
    ap.mass = 39.098;
    gaff_parameters.append(ap);
}
