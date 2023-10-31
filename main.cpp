#include "functions.h"
#include <string>
#include <fstream>
#include <complex>


void process(val_t eps, size_t n = 5)
{
    size_t k = 4;
    val_t h = (B-A)/k;
    val_t L = 1.6;


    for (int i = 1;;i++)
    {
        bool gauss = 1;
        val_t s1 = integrate(f,A,B,h,gauss,n);
        val_t s2 = integrate(f,A,B,h/L,gauss,n);
        val_t s3 = integrate(f,A,B,h/L/L,gauss,n);

        std::cout.precision(20);

        val_t m  = my_abs(my_log(my_abs((s3-s2)/(s2-s1)))/my_log(L));
        val_t cm = my_abs(s3-s2) / my_pow(h/L,m) / (1 - 1/my_pow(L,m));
        std::cout <<'#' << i << " Cm = " << cm << "  h = " << h/L/L << "  m = " << m << "  " << cm*my_pow(h/L/L, m)<<'\n'; 
        
        //std::cout.precision(std::numeric_limits<val_t>::digits10);
        std::cout.precision(20);
        std::cout << "Sh = " << s3 << "\n\n";
        std::cout.precision(3);

        if (cm*my_pow(h/L/L, m) < eps)
            break;

        //h /= 2;
        h = my_pow(eps/cm, 1/m) *L*L;
    }
}

signed main(int argc, char** argv)
{
    //std::cout.precision(std::numeric_limits<val_t>::digits10);
    process(1e-16, std::stoi(argv[1]));
    //std::cout << integrate_gauss(f,A,B,35).res << '\n';
    // 4.389308255579691622482841331293167979052130308324017142

    
    //
    //val_t exact("12.5818467721890266025251714551070231833486414946048158140881054801021790122047597920465776312160354"); 
    //val_t exact = std::stod("12.581846772189026602525");

    /*
    val_t minn = 1;
    val_t absmin;
    int min_n = -1;

    std::ofstream fout("data.txt");
    for (int n = 6; n < 50; n++)
    {
        auto [sum, abs_sum] = integrate_newton(f,A,B,n);
        fout << n << ' ' << my_log(abs_sum) << '\n';
        if (my_abs(sum - exact) < minn)
        {
            min_n = n;
            minn = my_abs(sum - exact);
            absmin = abs_sum;
        }
    }
    fout.close();
    std::cout << minn << '\n' << min_n << '\n' << absmin <<"\n\n";
    */
    //process(std::stod(argv[1]));
//    std::cout << integrate_gauss(f,A,B,std::stoi(argv[1])).res << '\n';
    //std::cout << integrate(f,A,B,hopt,true,m) << '\n';
    
}