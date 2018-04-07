/* Dates: 17/4/3-7 = Mon-Fri */
/* Add Execution Timing: 17/4/9 = Sun (around 13 sec) */
/* Author: Fu Lei */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define sq(x) ((x) * (x))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define FALSE 0
#define TRUE 1
#define DISCRETE 0
#define CONTINUOUS 1
#define EUROPEAN 0
#define AMERICAN 1
#define UP 1
#define MID 0
#define DOWN -1
#define QPRINT 0
#define RPRINT 1
#define BPRINT 2
#define VPRINT 3

typedef struct node Node, *pNode, **ppNode;

struct node {
    double q;
    double r;
    double b;
    double v;
    double pu;
    double pm;
    double pd;
    pNode u;
    pNode m;
    pNode d;
};

ppNode tree;

const double cr = 0.05,
             dT = 0.5,
             Pr = 100,
             kappa = 0.1,
             sigma = 0.015,
             T0 = 5,
             TN = 10,
             X = 100,
             y = 0.05;

int J, M, N;
double dt = 0.5, dr, *theta, *zero;

void mnset(void);
void zset(int);
void drset(void);
void jset(void);
void tset(void);
void rset(int);
void bset(void);
void vset(int);
double vget(void);
double vprice(double, int, int);
double prob(int, int, int);
double dotprod(double [], double [], int);
pNode node(int, int);
pNode u(int, int);
pNode m(int, int);
pNode d(int, int);
double *pu(int, int);
double *pm(int, int);
double *pd(int, int);
double *q(int, int);
double *r(int, int);
double *b(int, int);
double *v(int, int);
void aprint(double [], int);
void bprint(void);
void pprint(void);
void tprint(int);


int main(void) {

    clock_t tbegin = clock();

    printf("\nEuropean, Bounce (No-Bounce) if Interest Rate is Negative\n");
    for (double dtval = 0.5; dtval >= 1.0/500.0; dtval /= 2.0)
        printf("dt = %lg, N = %lg, value = %.4lf (%.4lf)\n",
                dtval, TN/dtval,
                vprice(dtval, TRUE, EUROPEAN),
                vprice(dtval, FALSE, EUROPEAN));

    printf("\nAmerican, Bounce (No-Bounce) if Interest Rate is Negative\n");
    for (double dtval = 0.5; dtval >= 1.0/500.0; dtval /= 2.0)
        printf("dt = %lg, N = %lg, value = %.4lf (%.4lf)\n",
                dtval, TN/dtval,
                vprice(dtval, TRUE, AMERICAN),
                vprice(dtval, FALSE, AMERICAN));

    clock_t tend = clock();

    printf("Execution time: %lf seconds\n", (double)(tend - tbegin) / CLOCKS_PER_SEC);

    return 0;
}

double vprice(double dtval, int bounce, int isAmerican) {
    dt = dtval;
    mnset();
    zset(DISCRETE);
    drset();
    jset();
    tset();
    rset(bounce);
    bset();
    vset(isAmerican);
    return vget();
}

void mnset(void) {
    M = T0 / dt;
    N = TN / dt;
}

void drset(void)
{   dr = sqrt(3*dt) * sigma; }

void zset(int ccompound) {
    zero = calloc(N+1, sizeof(double));
    for(int n = 0; n <= N; n++)
        zero[n] = ccompound ? exp(-1 * y * n * dt) : pow((1 + y * dT), -1 * (n * dt / dT));
}

void jset(void)
{   J = (int)(floor)(1 / (2 * kappa * dt)) + 1; }

void tset(void) {
    tree = (ppNode)calloc(N+1, sizeof(pNode));
    for (int n = 0; n <= N; n++) {
        int B = min(n,J);
        tree[n] = (pNode)calloc(2 * B + 1, sizeof(Node));
    }

    for (int n = 0; n < N; n++) {
        int B = min(n,J);
        for (int j = -1 * B; j <= B; j++) {
            int dj = (j == -1 * J) - (j == J);
            int k = j + dj;
            node(n,j)->u = node(n+1,k+1);
            node(n,j)->m = node(n+1,k  );
            node(n,j)->d = node(n+1,k-1);
            *pu(n,j) = prob(j, dj, UP);
            *pm(n,j) = prob(j, dj, MID);
            *pd(n,j) = prob(j, dj, DOWN);
        }
    }
}

void rset(int bounce) {
    theta = calloc(N-1, sizeof(double));
    *q(0,0) = 1;
    *r(0,0) = y;
    for (int n = 0; n < N - 1; n++) {
        int B = min(n,J);
        for (int j = -1 * B; j <= B; j++) {
            u(n,j)->q += *pu(n,j) * *q(n,j) * exp(-1 * *r(n,j) * dt);
            m(n,j)->q += *pm(n,j) * *q(n,j) * exp(-1 * *r(n,j) * dt);
            d(n,j)->q += *pd(n,j) * *q(n,j) * exp(-1 * *r(n,j) * dt);
        }
        B = min(n+1,J);
        double sum = 0;
        for (int j = -1 * B; j <= B; j++) {
            *r(n+1,j) = (1 - kappa * dt) * *r(n,0) + j * dr;
            sum += *q(n+1,j) * exp(-1 * *r(n+1,j) * dt);
        }
        theta[n] = log(sum / zero[n+2]) / (kappa * sq(dt));
        for (int j = -1 * B; j <= B; j++) {
            *r(n+1,j) += kappa * theta[n] * dt;
            if (bounce && *r(n+1,j) < 0) {
                node(n+1,j)->u = node(n+2,j+2);
                node(n+1,j)->m = node(n+2,j+1);
                node(n+1,j)->d = node(n+2,j);
                *pu(n+1,j) = prob(j, UP, UP);
                *pm(n+1,j) = prob(j, UP, MID);
                *pd(n+1,j) = prob(j, UP, DOWN);                
            }
        }
    }
}


void bset(void) {
    int B = min(N,J);
    for (int j = -1 * B; j <= B; j++)
        *b(N,j) = Pr * (1 + cr * dT);

    for (int n = N - 1; n >= 0; n--) {
        B = min(n,J);
        for (int j = -1 * B; j <= B; j++) {
            *b(n,j) = exp(-1 * *r(n,j) * dt) *
                      (*pu(n,j) * u(n,j)->b +
                       *pm(n,j) * m(n,j)->b +
                       *pd(n,j) * d(n,j)->b);
            if (fmod(n * dt, dT) == 0 && n > 0)
                *b(n,j) += Pr * cr * dT;
        }
    }
}

void vset(int isAmerican) {
    int B = min(M,J);

    for (int j = -1 * B; j <= B; j++) {
        *v(M,j) = *b(M,j) - Pr * cr * dT;
        *v(M,j) = *v(M,j) > X ? *v(M,j) - X : 0;
    }

    for (int n = M - 1; n >= 0; n--) {
        B = min(n,J);
        for (int j = -1 * B; j <= B; j++) {
            *v(n,j) = exp(-1 * *r(n,j) * dt) *
                      (*pu(n,j) * u(n,j)->v + 
                       *pm(n,j) * m(n,j)->v +
                       *pd(n,j) * d(n,j)->v);
        
            if (isAmerican) {
                double ap = fmod(n * dt, dT);
                if (ap == 0 && n > 0)
                    ap = dT;
                double ai = Pr * cr * ap;
                *v(n,j) = max(*v(n,j), *b(n,j) - ai - X);
            }
        }
    }
}

double vget(void)
{   return node(0,0)->v; }

double prob(int j, int dj, int move) {
    double a2 = j * kappa * dt + dj;
    double a1 = sq(a2);
    double a[3] = { 1, a1, a2 };
    double b[3];

    switch(move) {
        case UP   : memcpy(b, (double[]){  1.0/6.0,  0.5, -0.5 }, sizeof b); break;
        case MID  : memcpy(b, (double[]){  2.0/3.0, -1.0,  0.0 }, sizeof b); break;
        case DOWN : memcpy(b, (double[]){  1.0/6.0,  0.5,  0.5 }, sizeof b); break;
    }

    return dotprod(a, b, 3);
}

double dotprod(double a[], double b[], int dim) {
    double sum = 0.0;
    while (dim-- > 0)
        sum += a[dim] * b[dim];
    
    return sum;
}


/* The following functions are for data structure purposes */

pNode node(int n, int j)
{   return &tree[n][j + min(n,J)]; }

pNode u(int n, int j)
{   return tree[n][j + min(n,J)].u; }

pNode m(int n, int j)
{   return tree[n][j + min(n,J)].m; }

pNode d(int n, int j)
{   return tree[n][j + min(n,J)].d; }

double *pu(int n, int j)
{   return &tree[n][j+ min(n,J)].pu; }

double *pm(int n, int j)
{   return &tree[n][j+ min(n,J)].pm; }

double *pd(int n, int j)
{   return &tree[n][j+ min(n,J)].pd; }

double *q(int n, int j)
{   return &tree[n][j + min(n,J)].q; }

double *r(int n, int j)
{   return &tree[n][j + min(n,J)].r; }

double *b(int n, int j)
{   return &tree[n][j + min(n,J)].b; }

double *v(int n, int j)
{   return &tree[n][j + min(n,J)].v; }


/* The following functions are for debugging purposes only */

void aprint(double a[], int len) {
    for (int n = 0; n < len; n++) {
        printf("[%2d] %.2lf%%; ", n, a[n] * 100);
        if (n % 5 == 0)
            putchar('\n');
    }
}

void bprint(void) {
    for (int n = 0; n < N; n++) {
        int B = min(n,J);
        for (int j = B; j >= -1 * B; j--) {
            int dj = (j == -1 * J || *r(n,j) < 0) - (j == J);
            int k = j + dj;
            printf("(%2d,%3d) -> m(%2d,%3d) : %d, %d\n", n, j, n+1, k,
                   m(n,j) == node(n+1,k), *pm(n,j) == prob(j, dj, MID));
        }
        putchar('\n');
        
    }
}

void pprint(void) {
    for (int pattern = UP; pattern >= DOWN; pattern--)
        for (int j = -1 * J; j <= J; j++)
            for (int move = UP; move >= DOWN; move--) {
                char ptag, mtag;
                switch(pattern) {
                    case UP   : ptag = 'U'; break;
                    case MID  : ptag = 'M'; break;
                    case DOWN : ptag = 'D'; break;
                }
                switch(move) {
                    case UP   : mtag = 'u'; break;
                    case MID  : mtag = 'm'; break;
                    case DOWN : mtag = 'd'; break;
                }
                printf("%c: p%c(%d) = %8.2lf%%\n",ptag, mtag, j, prob(j, pattern, move) * 100);
            }
}

void tprint(int option) {
    int nPrint;
    switch(option) {
        case QPRINT : nPrint = N-1; break;
        case RPRINT : nPrint = N-1; break;
        case BPRINT : nPrint = N  ; break;
        case VPRINT : nPrint = M  ; break;
    }

    for (int n = 0; n <= nPrint; n++) {
        printf("\n[%d]\n", n);
        int B = min(n, J);
        for (int j = -1 * B; j <= B; j++) {
            double val;
            switch(option) {
                case QPRINT: val = *q(n,j); break;
                case RPRINT: val = *r(n,j); break;
                case BPRINT: val = *b(n,j); break;
                case VPRINT: val = *v(n,j); break;
            }
            char tag = j == 0 ? '*' : ' ';
            if (option == RPRINT)
	        printf("%8.2lf%c", val * 100, tag);
            else
                printf("%8.2lf%c", val, tag);
            if ((j + B) % 10 == 9 || j== B)
                putchar('\n');
        }
    }
}