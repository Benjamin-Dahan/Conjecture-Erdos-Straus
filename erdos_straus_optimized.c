#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

bool est_premier(int n) {
    if (n <= 1) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;

    int i;
    for (i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) {
            return false;
        }
    }
    return true;
}

bool check_conditions(int n) {
    if (!est_premier(n)) return false;
    if (n % 3 != 1 || n % 4 != 1 || n % 7 != 1 || (n % 7 != 1 && n % 7 != 4 && n % 7 != 2)) return false;

    int np3_div4 = (n + 3) / 4;
    if (np3_div4 % 2 == 0 || np3_div4 % 3 == 0 || np3_div4 % 5 == 0) return false;

    return true;
}

void main_loop() {
    int n, c, r, b;
    int sigma, csigma_minus_k;
    int stop = 0;  

    
    #pragma omp parallel for private(c, r, b, sigma, csigma_minus_k) schedule(dynamic)
    for (n = 3362; n < 10000000000; n++) {
        if (stop) continue;  

        if (check_conditions(n)) {
            for (c = 3; c < 200001; c++) {
                if (c % 4 == 3) {
                    r = n % c;

                    switch (r % 4) {
                        case 1: b = 0; break;
                        case 2: b = 1; break;
                        case 3: b = 2; break;
                        default: b = -1; continue;
                    }

                    sigma = (n + 4 * c - b * c - r) / (4 * c);
                    csigma_minus_k = c * sigma - (n + c) / 4;

                    if ((sigma * (n + c) / 4) % csigma_minus_k == 0) {
                        #pragma omp critical
                        {
                            printf("n=%d, (n+c)/4=%d, c=%d, sigma=%d\n", n, (n + c) / 4, c, sigma);
                        }
                        break; 
                    }
                }
            }
            if (c >= 200001) {
                #pragma omp critical
                {
                    printf("Le programme s'arrête : aucun c n'a fonctionné pour n=%d\n", n);
                }
                stop = 1; 
            }
        }
    }
}

int main() {
    main_loop();
    return 0;
}
