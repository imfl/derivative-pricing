/* Date: 17/4/9 = Sun */
/* Revised: 17/4/17 = Mon */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x) * (x))
#define TRUE 1
#define FALSE 0

typedef struct node Node, *pNode, **ppNode, ***pppNode;

struct node
{
    double *w;
    pNode uu;
    pNode du;
    pNode dd;
    pNode ud;
    pNode oo;
};

pppNode tree;

const double c = 0.04075, X = 0.87;

int N = 168, M = 21;

double dt, r, sigma1, sigma2, lambda, rho, u1, u2, X1, X2, prob[5];

void specify(void);
int build(void);
void terminal();
double backward();
double payoff(int, int, int, int);
int g(int, int, int, int);
double dotprod(double [], double [], int);
pNode node(int, int, int);
void parameters(void);

void specify(void)
{
    M = 10;
    N = M * 8;
    dt = 2.0 / N;
    r = 0.01;
    sigma1 = 0.05;
    sigma2 = 0.05;
    rho = -0.25;
    lambda = sqrt(5);
    u1 = exp(lambda * sigma1 * sqrt(dt));
    u2 = exp(lambda * sigma2 * sqrt(dt));
    X1 = log(X) / log(u1);
    X2 = log(X) / log(u2);
    double A = (sqrt(dt) / lambda) * (r - sq(sigma1)/2) / sigma1;
    double B = (sqrt(dt) / lambda) * (r - sq(sigma2)/2) / sigma2;
    double C = 1 / sq(lambda);
    double D = rho * C;
    double V [4] = {  A,  B,  C,  D };
    double U1[4] = {  1,  1,  1,  1 };
    double U2[4] = { -1,  1,  1, -1 };
    double U3[4] = { -1, -1,  1,  1 };
    double U4[4] = {  1, -1,  1, -1 };
    memcpy(prob, (double []) { dotprod(U1, V, 4) / 4,
                               dotprod(U2, V, 4) / 4,
                               dotprod(U3, V, 4) / 4,
                               dotprod(U4, V, 4) / 4,
                               1 - 1 / sq(lambda)  },
           sizeof prob);
}

int build(void)
{
    tree = (pppNode)calloc(N + 1, sizeof(ppNode));
    for (int n = 0; n <= N; n++) {
        tree[n] = (ppNode)calloc(2 * sq(n) + 2 * n + 1, sizeof(pNode));
        for (int p = 0; p < 2 * n + 1; p++)             
            tree[n][p] = (pNode)calloc(n + !(p % 2), sizeof(Node));
    }

    for (int n = 0; n <= N; n++)
        for (int i = -1 * n; i <= n; i++) {
            int J = n - (i + n) % 2; 
            for (int j = -1 * J; j <= J; j += 2) {
                pNode this = node(n,i,j);
                if (n < N) {
                    this->uu = node(n+1,i+1,j+1);
                    this->du = node(n+1,i-1,j+1);
                    this->dd = node(n+1,i-1,j-1);
                    this->ud = node(n+1,i+1,j-1);
                    this->oo = node(n+1,i  ,j  );
                }
                if (n <= M)
                    this->w = (double *)calloc(1, sizeof(double));
                else {
                    int m = n % M;
                    m = m ? m : M;
                    this->w = (double *)calloc(m + 1, sizeof(double));
                }
            }
        }

    return (N+1) * (2*sq(N) + 4*N + 3) / 3;    /* Number of nodes -- 3,217,929 in this case */
}

void terminal()
{
    for (int i = -1 * N; i <= N; i++) {
        int J = N - (i + N) % 2; 
        for (int j = -1 * J; j <= J; j += 2)
            for (int k = 0; k <= M; k++) {
                node(N,i,j)->w[k] = payoff(N, i, j, k);
            }
    }
    printf("Finish step %d\n", N);
}

double backward()
{
    for (int n = N - 1; n >= M; n--) {
        for (int i = -1 * n; i <= n; i++) {
            int J = n - (i + n) % 2; 
            for (int j = -1 * J; j <= J; j += 2) {
                pNode this = node(n,i,j);
                int m = n % M;
                m = (m || n == M) ? m : M;
                for (int k = 0; k <= m; k++) {
                    double next[5] = {
                        this->uu->w[g(n,i,j,k)],    /* At n = M (21), k shoule be M (21), because the   */
                        this->du->w[g(n,i,j,k)],    /* first period coupon is fixed. However, at that   */
                        this->dd->w[g(n,i,j,k)],    /* moment, the value of k does not affect grid      */
                        this->ud->w[g(n,i,j,k)],    /* function. So we pretend that k = 0 at that time. */
                        this->oo->w[g(n,i,j,k)]     /* By doing this we maintain consistency in form.   */
                    };
                    this->w[k] = exp(-1 * r * dt) * dotprod(next, prob, 5) +
                                 (n % M == 0) * c * (n == M ? M : k) / M;
                    if (i >= 0 && j >= 0 && (m == 0 || m == M))
                        this->w[k] = min(this->w[k], payoff(n,i,j,k));
                }
            }
        }
        printf("Finish step %d\n", n);
        for (int p = 0; p < 2 * n + 3; p++)
            free(tree[n+1][p]);
    }
    for (int n = M - 1; n >= 0; n--) {
        for (int i = -1 * n; i <= n; i++) {
            int J = n - (i + n) % 2;
            for (int j = -1 * J; j <= J; j += 2) {
                pNode this = node(n,i,j);
                double next[5] = {
                        this->uu->w[0],
                        this->du->w[0],
                        this->dd->w[0],
                        this->ud->w[0],
                        this->oo->w[0]
                    };
                this->w[0] = exp(-1 * r * dt) * (dotprod(next, prob, 5));
            }
        }
        printf("Finish step %d\n", n);
        for (int p = 0; p < 2 * n + 3; p++)
            free(tree[n+1][p]);
    }
    return node(0,0,0)->w[0];
}

double payoff(int n, int i, int j, int k)
{
    double h1 = i - X1, h2 = j - X2;

    return (double)c * k / M + 
           (h1 >= 0 && h2 >= 0) +
           (h1 < 0 || h2 < 0) * pow(h1 < h2 ? u1 : u2, min(h1, h2));
}

inline double dotprod(double a[], double b[], int I)
{
    double sum = 0;
    for (int i = 0; i < I; i++)
        sum += a[i] * b[i];

    return sum;
}

inline int g(int n, int i, int j, int k)
{
	return (i > X1 && j > X2) + k * (n % M != 0);
}

inline pNode node(int n, int i, int j)
{
	return &tree[n][i + n][(j + n - (i + n) % 2) / 2];
}

void parameters(void)
{
    printf("M = %d\n", M);
    printf("N = %d\n", N);
    printf("dt = %lf\n", dt);
    printf("r = %.2lf%%\n", r * 100);
    printf("sigma1 = %.2lf\n", sigma1);
    printf("sigma2 = %.2lf\n", sigma2);
    printf("rho = %.2lf\n", rho);
    printf("lambda = %.4lf\n", lambda);
    printf("u1 = %.4lf\n", u1);
    printf("u2 = %.4lf\n", u2);
    printf("X1 = %.4lf\n", X1);
    printf("X2 = %.4lf\n", X2);
    printf("puu = %.4lf\n", prob[0]);
    printf("pdu = %.4lf\n", prob[1]);
    printf("pdd = %.4lf\n", prob[2]);
    printf("pud = %.4lf\n", prob[3]);
    printf("poo = %.4lf\n", prob[4]);
}

void line()
{
    printf("-----------------------------\n");
}

int main(int argc, char *argv[])
{
	clock_t tbegin = clock();
	specify();
	build();
	terminal();
	double val = backward();
	clock_t tend = clock();
	specify();
    line();
	printf("Value of warrant = %.4lf\n", val);
	printf("Time elapsed = %.2lf seconds\n", (double)(tend - tbegin) / CLOCKS_PER_SEC);
	parameters();
    line();
	return 0;
}