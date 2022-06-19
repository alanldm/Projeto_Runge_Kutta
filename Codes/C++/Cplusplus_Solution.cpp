#include <iostream>
#include <list>
#include <math.h>

using namespace std;
double y_0 = -0.24593576;
double z_0 = -0.55769344;
double c[2] = {0, 0};
double t[2] = {0, 0};

double f(double x,double y){
  return -(100+(1.0/(x*x)))*y;
}

void coeficientes(double h, double x, double y, double z){
  double k1 = h*z;
  double l1 = h*f(x,y);
  double k2 = h*(z+l1/2.0);
  double l2 = h*f(x+h/2.0, y+k1/2.0);
  double k3 = h*(z+l2/2.0);
  double l3 = h*f(x+h/2.0, y+k2/2.0);
  double k4 = h*(z+l3);
  double l4 = h*f(x+h, y+k3);

  c[0] = l1+2*l2+2*l3+l4;
  c[1] = k1+2*k2+2*k3+k4;
}

void termos(double h, double x, double y, double z){
  double k0 = (pow(h,2))*f(x, y);
  double k1 = (pow(h,2))*f(x+(h/4.0), y+(h*z/4.0)+(k0/32.0));
  double k2 = (pow(h,2))*f(x+(h/2.0), y+(h*z/2.0)-(k0/24.0)+(k1/6.0));
  double k3 = (pow(h,2))*f(x+(3.0*h/4), y+(3.0*h*z/4)+(3.0*k0/32)+(k1/8.0)+(k2/16.0));
  double k4 = (pow(h,2))*f(x+h, y+h*z+(3.0*k1/7)-(k2/14.0)+(k3/7.0));

  t[0] = 7*k0+32*k1+12*k2+32*k3+7*k4;
  t[1] = 7*k0+24*k1+6*k2+8*k3;
}

int main() {

  double h[3] = {0.025, 0.25, 0.5};

  for(int i=0; i<3; i++){
    if(h[i]==0.025){
      double x_025[1217];
      double y_rk4_025[1217];
      double z_rk4_025[1217];
      double y_rk6_025[1217];
      double z_rk6_025[1217];
      
      y_rk4_025[0] = y_0;
      z_rk4_025[0] = z_0;
      y_rk6_025[0] = y_0;
      z_rk6_025[0] = z_0; 
        
      double cont = 1;
      for(int j=0; j<1217; j++){
        x_025[j] = cont;
        cont += 0.025;
      }

      for(int n=1; n<1217; n++){
        coeficientes(0.025, x_025[n-1], y_rk4_025[n-1], z_rk4_025[n-1]);
        z_rk4_025[n] = z_rk4_025[n-1] + (1.0/6)*c[0];
        y_rk4_025[n] = y_rk4_025[n-1] + (1.0/6)*c[1];

        termos(0.025, x_025[n-1], y_rk6_025[n-1], z_rk6_025[n-1]);
        z_rk6_025[n] = z_rk6_025[n-1] + (1.0/(90*0.025))*t[0];
        y_rk6_025[n] = y_rk6_025[n-1] + 0.025*z_rk6_025[n-1] + (1.0/90)*t[1];
      }
      cout << "RK4" << "     0.025     " << "RK6" << endl;
      for(int n=0; n<1217; n++){
        cout << y_rk4_025[n] << "          " << y_rk6_025[n] << endl;
      }
    }else if(h[i]==0.25){
      double x_25[122];
      double y_rk4_25[122];
      double z_rk4_25[122];
      double y_rk6_25[122];
      double z_rk6_25[122];
      
      y_rk4_25[0] = y_0;
      z_rk4_25[0] = z_0;
      y_rk6_25[0] = y_0;
      z_rk6_25[0] = z_0; 
        
      double cont = 1;
      for(int j=0; j<122; j++){
        x_25[j] = cont;
        cont += 0.25;
      }

      for(int n=1; n<122; n++){
        coeficientes(0.25, x_25[n-1], y_rk4_25[n-1], z_rk4_25[n-1]);
        z_rk4_25[n] = z_rk4_25[n-1] + (1.0/6)*c[0];
        y_rk4_25[n] = y_rk4_25[n-1] + (1.0/6)*c[1];

        termos(0.25, x_25[n-1], y_rk6_25[n-1], z_rk6_25[n-1]);
        z_rk6_25[n] = z_rk6_25[n-1] + (1.0/(90*0.25))*t[0];
        y_rk6_25[n] = y_rk6_25[n-1] + 0.25*z_rk6_25[n-1] + (1.0/90)*t[1];
      }
      cout << "RK4" << "     0.25     " << "RK6" << endl;
      for(int n=0; n<122; n++){
        cout << y_rk4_25[n] << "          " << y_rk6_25[n] << endl;
      }
    }else{
      double x_5[61];
      double y_rk4_5[61];
      double z_rk4_5[61];
      double y_rk6_5[61];
      double z_rk6_5[61];

      for(int n=0; n<61; n++){
        y_rk4_5[n] = 0;
        z_rk4_5[n] = 0;
        y_rk6_5[n] = 0;
        z_rk6_5[n] = 0; 
      }
      y_rk4_5[0] = y_0;
      z_rk4_5[0] = z_0;
      y_rk6_5[0] = y_0;
      z_rk6_5[0] = z_0; 
        
      double cont = 1;
      for(int j=0; j<61; j++){
        x_5[j] = cont;
        cont += 0.5;
      }

      for(int n=1; n<61; n++){
        coeficientes(0.5, x_5[n-1], y_rk4_5[n-1], z_rk4_5[n-1]);
        z_rk4_5[n] = z_rk4_5[n-1] + (1.0/6)*c[0];
        y_rk4_5[n] = y_rk4_5[n-1] + (1.0/6)*c[1];        

        termos(0.5, x_5[n-1], y_rk6_5[n-1], z_rk6_5[n-1]);
        z_rk6_5[n] = z_rk6_5[n-1] + (1.0/(90*0.5))*t[0];
        y_rk6_5[n] = y_rk6_5[n-1] + 0.5*z_rk6_5[n-1] + (1.0/90)*t[1];
      }
      cout << "RK4" << "     0.5     " << "RK6" << endl;
      for(int n=0; n<61; n++){
        cout << y_rk4_5[n] << "          " << y_rk6_5[n] << endl;
      }
    }
  }
}