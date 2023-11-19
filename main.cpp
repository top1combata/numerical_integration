#include "functions.h"
#include <string>
#include <fstream>


signed main(int argc, char** argv)
{
    using namespace std;
    std::function<val_t(val_t)> f = [](val_t x)->val_t {return 4*cos(x/2)*exp(-2*x/3) + val_t(12)/5 * sin(9*x/2)*exp(x/8)+2;};
    //std::function<val_t(val_t)> f = [](val_t x)->val_t {return pow(x,10)+pow(x,9)+pow(x,8)+x*x+x+1+exp(x);};
    val_t A = val_t(13)/10, B = val_t(22)/10, beta = val_t(5)/6;

    
    //std::cout.precision(std::numeric_limits<val_t>::digits10);
    string s("12.5818467721890266025251714551070231833486414946048158140881054801021790122047597920465776312160354");
    //string s("23.17747840368038699937255372663352073488484164038818032697609060901458122750935151781353633120826208");
#ifdef MULTIPRECISION
    val_t exact(s);
#else
    val_t exact = stold(s);
#endif

    int k = 3;
    val_t h = 0.1;
    val_t m,cm;
    val_t L = 1.5;
    val_t sum1, sum2 = integrate_qaws(f,A,B,0,beta,h,1,k);
    h /= L;
    val_t sum3 = integrate_qaws(f,A,B,0,beta,h,1,k);

    printf("%6s %7s %13s %10s %10s\n", "h", "m", "cm", "est err","abs err");

    for (int i = 0; ; i++)
    {
        h /= L;
        sum1 = sum2;
        sum2 = sum3;
        sum3 = integrate_qaws(f,A,B,0,beta,h,1,k);
        

        m = abs(log(abs(sum3-sum2)/abs(sum2-sum1)) / log(L));
        cm = abs(sum3-sum2)/pow(h,m)/(pow(L,m)-1);

        printf("%10.3e %7.3f %10.3e %10.3e %10.3e\n", double(h), double(m), double(cm), double(cm*pow(h,m)), double(abs(exact-sum3)));
    }
    
    return 0;
}