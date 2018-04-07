/* Date: 17/4/9 = Sun */
/* Revise: 17/4/17 = Mon */
/* Review: 17/4/24 = Mon */

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

struct node {
    double *w;
    pNode uu;
    pNode du;
    pNode dd;
    pNode ud;
    pNode oo;
};

pppNode tree;

const double c = 0.04075, X = 0.87;

int N = 168, M = 21, callable = TRUE;

double dt, r, sigma1, sigma2, lambda, rho, u1, u2, x1, x2, prob[5];

void compute(int, double, double, double, double, int);
void value(int, double, double, double, double, int);
void specify(int, double, double, double, double, int);
int check(void);
int build(void);
void terminal(void);
double backward(void);
double payoff(int, int, int, int);
int g(int, int, int, int);
double dotprod(double [], double [], int);
pNode node(int, int, int);
void parameters(int, int, int, int);

void compute(int v_M, double v_r, double v_sigma1, double v_sigma2, double v_rho, int v_callable) {
    
    specify(v_M, v_r, v_sigma1, v_sigma2, v_rho, v_callable);
    
    parameters(TRUE, TRUE, TRUE, TRUE);
    
    printf("\n");
}

void value(int v_M, double v_r, double v_sigma1, double v_sigma2, double v_rho, int v_callable) {
    
    specify(v_M, v_r, v_sigma1, v_sigma2, v_rho, v_callable);

    if (check())
        return;

    clock_t tbegin = clock();

    build();

    terminal();

    double val = backward();

    clock_t tend = clock();

    parameters(TRUE, TRUE, TRUE, FALSE);

    printf("value of warrnt = %8.6lf | ", val);
   
    printf("time elapsed = %.2lf seconds\n", (double)(tend - tbegin) / CLOCKS_PER_SEC);
}


void specify(int v_M, double v_r, double v_sigma1, double v_sigma2, double v_rho, int v_callable) {

    M = v_M;
    N = M * 8;
    dt = 2.0 / N;
    r = v_r;
    sigma1 = v_sigma1;
    sigma2 = v_sigma2;
    rho = v_rho;
    callable = v_callable;
    lambda = sqrt(3);

    u1 = exp(lambda * sigma1 * sqrt(dt));
    u2 = exp(lambda * sigma2 * sqrt(dt));
    x1 = log(X) / log(u1);
    x2 = log(X) / log(u2);

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

int check(void) {

    for (int i = 0; i < 5; i++)
        if (prob[i] < 0) {
            parameters(TRUE, TRUE, TRUE, TRUE);
            printf("\n");
            return TRUE;
        }

    return FALSE;
}

int build(void) {
    
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

void terminal(void) {
    
    for (int i = -1 * N; i <= N; i++) {
        int J = N - (i + N) % 2; 
        for (int j = -1 * J; j <= J; j += 2)
            for (int k = 0; k <= M; k++)
                node(N,i,j)->w[k] = payoff(N, i, j, k);
    }
}

double backward(void) {

    for (int n = N - 1; n >= M; n--) {
        for (int i = -1 * n; i <= n; i++) {
            int J = n - (i + n) % 2; 
            for (int j = -1 * J; j <= J; j += 2) {
                pNode this = node(n,i,j);
                int m = n % M;
                m = (m || n == M) ? m : M;
                for (int k = 0; k <= m; k++) {
                    double next[5] = {
                        this->uu->w[g(n,i,j,k)],   /* At n = M (21), k shoule be M (21), because the   */
                        this->du->w[g(n,i,j,k)],   /* first period coupon is fixed. However, at that   */
                        this->dd->w[g(n,i,j,k)],   /* moment, the value of k does not affect grid      */
                        this->ud->w[g(n,i,j,k)],   /* function. So we pretend that k = 0 at that time. */
                        this->oo->w[g(n,i,j,k)]    /* By doing this we maintain consistency in form.   */
                    };
                    this->w[k] = exp(-1 * r * dt) * dotprod(next, prob, 5) +
                                 (n % M == 0) * c * (n == M ? M : k) / M;
                    if (callable && i >= 0 && j >= 0 && (m == 0 || m == M))
                        this->w[k] = min(this->w[k], payoff(n,i,j,k));
                }
            }
        }
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
                this->w[0] = exp(-1 * r * dt) * dotprod(next, prob, 5);
            }
        }
        for (int p = 0; p < 2 * n + 3; p++)
            free(tree[n+1][p]);
    }

    double warrant = node(0,0,0)->w[0];

    free(tree[0][0]);

    free(tree[0]);

    free(tree);

    return warrant;
}

double payoff(int n, int i, int j, int k) {
    
    double h1 = i - x1, h2 = j - x2;

    return (double)c * k / M + 
           (h1 >= 0 && h2 >= 0) +
           (h1 < 0 || h2 < 0) * ((sigma1 * h1 < sigma2 * h2) ? pow(u1, h1) : pow(u2, h2));
}

inline double dotprod(double a[], double b[], int I) {
    
    double sum = 0;
    
    for (int i = 0; i < I; i++)
        sum += a[i] * b[i];

    return sum;
}

inline int g(int n, int i, int j, int k)
{   return (i > x1 && j > x2) + k * (n % M != 0); }

inline pNode node(int n, int i, int j)
{   return &tree[n][i + n][(j + n - (i + n) % 2) / 2]; }

void parameters(int p_r, int p_sigma, int p_rho, int p_prob) {

    if (p_r)
        printf("r = %4.1lf%% | ", r * 100);

    if (p_sigma) {
        printf("sigma1 = %4.2lf | ", sigma1);
        printf("sigma2 = %4.2lf | ", sigma2);
    }

    if (p_rho)
        printf("rho = %5.2lf | ", rho);

    if (p_prob) {
        printf("puu %c %5.2lf | ", prob[0] < 0 ? '*' : '=', prob[0]);
        printf("pdu %c %5.2lf | ", prob[1] < 0 ? '*' : '=', prob[1]);
        printf("pdd %c %5.2lf | ", prob[2] < 0 ? '*' : '=', prob[2]);
        printf("pud %c %5.2lf | ", prob[3] < 0 ? '*' : '=', prob[3]);
        printf("poo %c %5.2lf",    prob[4] < 0 ? '*' : '=', prob[4]);
    }
}

int main(int argc, char *argv[]) {

    printf("Enter 0 to value warrants, 1 to compute probabilities: ");

    int flag = TRUE;
    scanf("%d", &flag);
    void (* proceed)(int, double, double, double, double, int);
    proceed = flag ? compute : value;

    int    v_M = 21,
           v_callable = TRUE;
    double      r_from,	     r_to,   	r_by,
           sigma1_from, sigma1_to, sigma1_by,
           sigma2_from, sigma2_to, sigma2_by,
	      rho_from,	   rho_to,    rho_by;

    printf("M = ");					scanf( "%d", &v_M);
    printf("r ranges from = ");				scanf("%lf", &r_from);
    printf("r ranges to = ");				scanf("%lf", &r_to);
    printf("r ranges by = ");				scanf("%lf", &r_by);
    printf("sigma1 ranges from = ");			scanf("%lf", &sigma1_from);
    printf("sigma1 ranges to = ");			scanf("%lf", &sigma1_to);
    printf("sigma1 ranges by = ");			scanf("%lf", &sigma1_by);
    printf("sigma2 ranges from = ");			scanf("%lf", &sigma2_from);
    printf("sigma2 ranges to = ");			scanf("%lf", &sigma2_to);
    printf("sigma2 ranges by = ");			scanf("%lf", &sigma2_by);
    printf("rho ranges from = ");			scanf("%lf", &rho_from);
    printf("rho ranges to = ");				scanf("%lf", &rho_to);
    printf("rho ranges by = ");				scanf("%lf", &rho_by);
    printf("callable (1 for true, 0 for false) = ");	scanf( "%d", &v_callable);

    for (double v_r 	 = 	r_from; v_r	 <= 	 r_to; v_r	+= 	r_by)
    for (double v_sigma1 = sigma1_from; v_sigma1 <= sigma1_to; v_sigma1 += sigma1_by)
    for (double v_sigma2 = sigma2_from; v_sigma2 <= sigma2_to; v_sigma2 += sigma2_by)
    for (double v_rho 	 =    rho_from; v_rho	 <=    rho_to; v_rho	+=    rho_by)
    	
        (*proceed)(v_M, v_r, v_sigma1, v_sigma2, v_rho, v_callable);

    return 0;
}
